#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <vulkan/vulkan.h>

namespace bc {

// Simple BC1 block encoder
inline void encodeBC1Block(const uint8_t* rgba, uint8_t* block, int stride) {
    // Find min/max colors in the 4x4 block
    uint8_t minR = 255, minG = 255, minB = 255;
    uint8_t maxR = 0, maxG = 0, maxB = 0;

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint8_t* pixel = rgba + (y * stride + x * 4);
            minR = std::min(minR, pixel[0]);
            minG = std::min(minG, pixel[1]);
            minB = std::min(minB, pixel[2]);
            maxR = std::max(maxR, pixel[0]);
            maxG = std::max(maxG, pixel[1]);
            maxB = std::max(maxB, pixel[2]);
        }
    }

    // Convert to RGB565
    uint16_t color0 = ((maxR >> 3) << 11) | ((maxG >> 2) << 5) | (maxB >> 3);
    uint16_t color1 = ((minR >> 3) << 11) | ((minG >> 2) << 5) | (minB >> 3);

    // Ensure color0 > color1 for opaque mode
    if (color0 < color1) {
        std::swap(color0, color1);
        std::swap(minR, maxR);
        std::swap(minG, maxG);
        std::swap(minB, maxB);
    }

    // Write colors
    block[0] = color0 & 0xFF;
    block[1] = (color0 >> 8) & 0xFF;
    block[2] = color1 & 0xFF;
    block[3] = (color1 >> 8) & 0xFF;

    // Calculate interpolated colors
    uint8_t colors[4][3];
    colors[0][0] = maxR; colors[0][1] = maxG; colors[0][2] = maxB;
    colors[1][0] = minR; colors[1][1] = minG; colors[1][2] = minB;
    colors[2][0] = (2 * maxR + minR) / 3;
    colors[2][1] = (2 * maxG + minG) / 3;
    colors[2][2] = (2 * maxB + minB) / 3;
    colors[3][0] = (maxR + 2 * minR) / 3;
    colors[3][1] = (maxG + 2 * minG) / 3;
    colors[3][2] = (maxB + 2 * minB) / 3;

    // Choose best color for each pixel
    uint32_t indices = 0;
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint8_t* pixel = rgba + (y * stride + x * 4);

            int bestIdx = 0;
            int bestDist = INT32_MAX;

            for (int i = 0; i < 4; i++) {
                int dr = (int)pixel[0] - (int)colors[i][0];
                int dg = (int)pixel[1] - (int)colors[i][1];
                int db = (int)pixel[2] - (int)colors[i][2];
                int dist = dr*dr + dg*dg + db*db;

                if (dist < bestDist) {
                    bestDist = dist;
                    bestIdx = i;
                }
            }

            indices |= (bestIdx << ((y * 4 + x) * 2));
        }
    }

    // Write indices
    block[4] = indices & 0xFF;
    block[5] = (indices >> 8) & 0xFF;
    block[6] = (indices >> 16) & 0xFF;
    block[7] = (indices >> 24) & 0xFF;
}

// Simple BC2 block encoder (explicit 4-bit alpha + BC1 color)
inline void encodeBC2Block(const uint8_t* rgba, uint8_t* block, int stride) {
    // Alpha block (first 8 bytes) - 4 bits per pixel
    for (int y = 0; y < 4; y++) {
        uint16_t row = 0;
        for (int x = 0; x < 4; x++) {
            const uint8_t* pixel = rgba + (y * stride + x * 4);
            uint8_t alpha = pixel[3];
            // Convert 8-bit alpha to 4-bit
            uint8_t alpha4 = alpha >> 4;
            row |= (alpha4 << (x * 4));
        }
        block[y * 2] = row & 0xFF;
        block[y * 2 + 1] = (row >> 8) & 0xFF;
    }

    // Color block (last 8 bytes) - use BC1 encoding
    encodeBC1Block(rgba, block + 8, stride);
}

// Simple BC3 block encoder (interpolated alpha + BC1 color)
inline void encodeBC3Block(const uint8_t* rgba, uint8_t* block, int stride) {
    // Alpha block (first 8 bytes)
    uint8_t minA = 255, maxA = 0;

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint8_t* pixel = rgba + (y * stride + x * 4);
            minA = std::min(minA, pixel[3]);
            maxA = std::max(maxA, pixel[3]);
        }
    }

    block[0] = maxA;
    block[1] = minA;

    // Calculate alpha indices
    uint64_t alphaIndices = 0;
    int bitPos = 0;

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint8_t* pixel = rgba + (y * stride + x * 4);
            uint8_t alpha = pixel[3];

            int idx;
            if (maxA == minA) {
                idx = 0;
            } else {
                // Simple linear interpolation
                int range = maxA - minA;
                idx = ((alpha - minA) * 7 + range/2) / range;
                idx = std::min(7, std::max(0, idx));
            }

            alphaIndices |= (uint64_t(idx) << bitPos);
            bitPos += 3;
        }
    }

    // Write alpha indices (48 bits)
    for (int i = 0; i < 6; i++) {
        block[2 + i] = (alphaIndices >> (i * 8)) & 0xFF;
    }

    // Color block (last 8 bytes) - use BC1 encoding
    encodeBC1Block(rgba, block + 8, stride);
}

// Main compression function
inline std::vector<uint8_t> compressToBC(VkFormat format, const uint8_t* rgbaData,
                                         uint32_t width, uint32_t height) {
    uint32_t blocksX = (width + 3) / 4;
    uint32_t blocksY = (height + 3) / 4;

    size_t blockSize = 8;  // BC1 and BC4
    bool useBC3 = false;

    switch (format) {
        case VK_FORMAT_BC1_RGB_UNORM_BLOCK:
        case VK_FORMAT_BC1_RGB_SRGB_BLOCK:
        case VK_FORMAT_BC1_RGBA_UNORM_BLOCK:
        case VK_FORMAT_BC1_RGBA_SRGB_BLOCK:
            blockSize = 8;
            break;

        case VK_FORMAT_BC2_UNORM_BLOCK:
        case VK_FORMAT_BC2_SRGB_BLOCK:
            blockSize = 16;
            // BC2 will use its own encoder
            break;

        case VK_FORMAT_BC3_UNORM_BLOCK:
        case VK_FORMAT_BC3_SRGB_BLOCK:
            blockSize = 16;
            useBC3 = true;
            break;

        case VK_FORMAT_BC4_UNORM_BLOCK:
        case VK_FORMAT_BC4_SNORM_BLOCK:
            blockSize = 8;
            break;

        case VK_FORMAT_BC5_UNORM_BLOCK:
        case VK_FORMAT_BC5_SNORM_BLOCK:
        case VK_FORMAT_BC6H_UFLOAT_BLOCK:
        case VK_FORMAT_BC6H_SFLOAT_BLOCK:
        case VK_FORMAT_BC7_UNORM_BLOCK:
        case VK_FORMAT_BC7_SRGB_BLOCK:
            blockSize = 16;
            break;

        default:
            return {};  // Unsupported format
    }

    std::vector<uint8_t> compressed(blocksX * blocksY * blockSize);

    for (uint32_t by = 0; by < blocksY; by++) {
        for (uint32_t bx = 0; bx < blocksX; bx++) {
            uint8_t* blockDst = compressed.data() + (by * blocksX + bx) * blockSize;

            // Create a 4x4 RGBA block with clamping at edges
            uint8_t blockRGBA[64];  // 4x4x4 bytes
            for (int py = 0; py < 4; py++) {
                for (int px = 0; px < 4; px++) {
                    uint32_t sy = std::min(by * 4 + py, height - 1);
                    uint32_t sx = std::min(bx * 4 + px, width - 1);
                    const uint8_t* src = rgbaData + (sy * width + sx) * 4;
                    uint8_t* dst = blockRGBA + (py * 4 + px) * 4;
                    dst[0] = src[0];
                    dst[1] = src[1];
                    dst[2] = src[2];
                    dst[3] = src[3];
                }
            }

            // Encode the block
            if (format == VK_FORMAT_BC2_UNORM_BLOCK || format == VK_FORMAT_BC2_SRGB_BLOCK) {
                encodeBC2Block(blockRGBA, blockDst, 16);  // 4 pixels * 4 bytes = 16 byte stride
            } else if (useBC3) {
                encodeBC3Block(blockRGBA, blockDst, 16);
            } else {
                encodeBC1Block(blockRGBA, blockDst, 16);
            }
        }
    }

    return compressed;
}

} // namespace bc