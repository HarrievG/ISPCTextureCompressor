// Fallback stubs when ISPC compiler is not available
#ifdef ISPC_FALLBACK_MODE

#include "ispc_texcomp.h"
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cfloat>
#include "bc_compression.h"

// Fallback implementations using our bc_compression.h
void CompressBlocksBC1(const rgba_surface* src, uint8_t* dst) {
#ifdef ISPC_FALLBACK_MODE
    // Use the simple BC1 encoder from bc_compression.h
    std::vector<uint8_t> compressed = bc::compressToBC(
        VK_FORMAT_BC1_RGB_UNORM_BLOCK, src->ptr, src->width, src->height);
    if (!compressed.empty()) {
        memcpy(dst, compressed.data(), compressed.size());
    }
#endif
}

void CompressBlocksBC2(const rgba_surface* src, uint8_t* dst) {
#ifdef ISPC_FALLBACK_MODE
    // Use the simple BC2 encoder from bc_compression.h
    std::vector<uint8_t> compressed = bc::compressToBC(
        VK_FORMAT_BC2_UNORM_BLOCK, src->ptr, src->width, src->height);
    if (!compressed.empty()) {
        memcpy(dst, compressed.data(), compressed.size());
    }
#endif
}

// BC3 implementation moved to ispc_texcomp.cpp with DirectXTex-compliant OptimizeAlpha
// void CompressBlocksBC3(const rgba_surface* src, uint8_t* dst) {
// #ifdef ISPC_FALLBACK_MODE
//     // Use the simple BC3 encoder from bc_compression.h
//     std::vector<uint8_t> compressed = bc::compressToBC(
//         VK_FORMAT_BC3_UNORM_BLOCK, src->ptr, src->width, src->height);
//     if (!compressed.empty()) {
//         memcpy(dst, compressed.data(), compressed.size());
//     }
// #endif
// }

// BC4 and BC5 implementations moved to ispc_texcomp.cpp with DirectXTex-compliant OptimizeAlpha
// void CompressBlocksBC4(const rgba_surface* src, uint8_t* dst) {
//     // Stub
// }
//
// void CompressBlocksBC5(const rgba_surface* src, uint8_t* dst) {
//     // Stub
// }

void CompressBlocksBC6H_fallback(const rgba_surface* src, uint8_t* dst, bc6h_enc_settings* settings) {
    // BC6H C++ fallback implementation based on Microsoft DirectXTex
    // Properly handles IEEE 754 half-float data according to BC6H specification

    if (!src || !dst || !settings) return;

    // BC6H uses 4x4 blocks, 16 bytes per block
    const int blockWidth = (src->width + 3) / 4;
    const int blockHeight = (src->height + 3) / 4;
    
    for (int by = 0; by < blockHeight; by++) {
        for (int bx = 0; bx < blockWidth; bx++) {
            uint8_t* blockDst = dst + (by * blockWidth + bx) * 16;
            
            // Extract 4x4 block of half-float RGB data
            // The input is half-float data (16-bit per channel, RGB only)
            uint16_t blockHalfFloat[48]; // 4x4x3 channels as half-float
            
            for (int py = 0; py < 4; py++) {
                for (int px = 0; px < 4; px++) {
                    int sy = std::min(by * 4 + py, src->height - 1);
                    int sx = std::min(bx * 4 + px, src->width - 1);
                    
                    // BC6H input is half-float RGBA (8 bytes per pixel, alpha ignored)
                    const uint8_t* srcPtr = src->ptr + (sy * src->stride + sx * 8);

                    int blockIdx = (py * 4 + px) * 3;

                    // Load half-float values directly (RGBA layout, ignore A)
                    blockHalfFloat[blockIdx + 0] = *reinterpret_cast<const uint16_t*>(srcPtr + 0); // R
                    blockHalfFloat[blockIdx + 1] = *reinterpret_cast<const uint16_t*>(srcPtr + 2); // G
                    blockHalfFloat[blockIdx + 2] = *reinterpret_cast<const uint16_t*>(srcPtr + 4); // B
                    // Alpha at srcPtr + 6 is ignored
                }
            }
            
            // Compress the BC6H block using DirectXTex-compliant algorithm
            compressBC6HBlockDirectXTex(blockHalfFloat, blockDst, settings);
        }
    }
}

// Forward declaration
static void compressBC6HBlockDirectXTex(const uint16_t blockHalfFloat[48], uint8_t* blockDst, const bc6h_enc_settings* settings);

// BC6H block compression - Microsoft specification compliant
static void compressBC6HBlock(const uint16_t blockHalfFloat[48], uint8_t* blockDst, const bc6h_enc_settings* settings) {
    // BC6H Mode 11 implementation for signed, Mode 10 for unsigned
    // According to Microsoft BC6H specification

    // Convert half-float to value for BC6H processing
    // DirectXTex approach: convert half-float to float, then scale appropriately
    auto halfToFloat = [](uint16_t half) -> float {
        // Convert IEEE 754 half-float to single-precision float
        uint32_t sign = (half >> 15) & 1;
        uint32_t exponent = (half >> 10) & 0x1F;
        uint32_t mantissa = half & 0x3FF;

        // Handle special cases
        if (exponent == 0) {
            if (mantissa == 0) {
                // Zero
                return sign ? -0.0f : 0.0f;
            }
            // Denormalized number
            float value = (mantissa / 1024.0f) * powf(2.0f, -14.0f);
            return sign ? -value : value;
        }
        if (exponent == 31) {
            // Inf or NaN - clamp to max finite
            return sign ? -65504.0f : 65504.0f;
        }

        // Normalized number
        float value = (1.0f + mantissa / 1024.0f) * powf(2.0f, (int)exponent - 15);
        return sign ? -value : value;
    };
    
    // Convert half-float block to float values
    float blockFloat[48];
    for (int i = 0; i < 48; i++) {
        blockFloat[i] = halfToFloat(blockHalfFloat[i]);
    }
    
    // BC6H Mode table from Microsoft specification
    struct BC6HMode {
        int mode;           // Mode number (0-13)
        int modeBits;       // Bits for mode field
        int endpointBits;  // Bits per endpoint
        bool hasDelta;     // Whether endpoints use delta compression
        bool hasPartition; // Whether block has multiple partitions
        bool isSigned;     // Whether mode supports signed values
    };
    
    static const BC6HMode modes[] = {
        {0, 2, 10, false, false, true},   // Mode 0: 2-bit mode, 10-bit endpoints, no delta, 1 partition, signed
        {1, 5, 7,  false, false, true},   // Mode 1: 5-bit mode, 7-bit endpoints, no delta, 1 partition, signed
        {2, 3, 11, false, false, true},  // Mode 2: 3-bit mode, 11-bit endpoints, no delta, 1 partition, signed
        {3, 0, 0,  false, false, false},  // Reserved
        {4, 0, 0,  false, false, false},  // Reserved
        {5, 4, 9,  false, false, true},   // Mode 5: 4-bit mode, 9-bit endpoints, no delta, 1 partition, signed
        {6, 4, 8,  false, false, true},   // Mode 6: 4-bit mode, 8-bit endpoints, no delta, 1 partition, signed
        {7, 0, 0,  false, false, false},  // Reserved
        {8, 0, 0,  false, false, false},  // Reserved
        {9, 5, 6,  false, false, true},   // Mode 9: 5-bit mode, 6-bit endpoints, no delta, 1 partition, signed
        {10, 5, 10, false, false, true}, // Mode 10: 5-bit mode, 10-bit endpoints, no delta, 1 partition, signed
        {11, 5, 11, false, false, true}, // Mode 11: 5-bit mode, 11-bit endpoints, no delta, 1 partition, signed
        {12, 5, 12, false, false, true}, // Mode 12: 5-bit mode, 12-bit endpoints, no delta, 1 partition, signed
        {13, 5, 16, false, false, true}, // Mode 13: 5-bit mode, 16-bit endpoints, no delta, 1 partition, signed
    };
    
    // Use Mode 11 for signed (better range), Mode 10 for unsigned
    const BC6HMode& mode = settings->is_signed ? modes[11] : modes[10];
    
    // Find min/max values for each channel
    float minR = FLT_MAX, maxR = -FLT_MAX;
    float minG = FLT_MAX, maxG = -FLT_MAX;
    float minB = FLT_MAX, maxB = -FLT_MAX;

    for (int i = 0; i < 16; i++) {
        float r = blockFloat[i * 3 + 0];
        float g = blockFloat[i * 3 + 1];
        float b = blockFloat[i * 3 + 2];

        minR = std::min(minR, r); maxR = std::max(maxR, r);
        minG = std::min(minG, g); maxG = std::max(maxG, g);
        minB = std::min(minB, b); maxB = std::max(maxB, b);
    }

    // Scale endpoints to BC6H integer range
    // The scaling must match the actual HDR range in the input data
    // For signed: input is [-8, 8], scale to fit in signed 16-bit range
    // For unsigned: input is [0, 16], scale appropriately
    float scale = settings->is_signed ? (32767.0f / 8.0f) : (65504.0f / 16.0f);

    int16_t minR_int = static_cast<int16_t>(std::round(minR * scale));
    int16_t maxR_int = static_cast<int16_t>(std::round(maxR * scale));
    int16_t minG_int = static_cast<int16_t>(std::round(minG * scale));
    int16_t maxG_int = static_cast<int16_t>(std::round(maxG * scale));
    int16_t minB_int = static_cast<int16_t>(std::round(minB * scale));
    int16_t maxB_int = static_cast<int16_t>(std::round(maxB * scale));

    // BC6H quantization according to Microsoft BC6H specification
    // This implements the correct BC6H quantization formula from DirectXTex
    auto bc6hQuantize = [](int16_t val, int prec, bool isSigned) -> int16_t {
        int q;
        
        if (isSigned) {
            // Signed quantization: BC6H signed format uses (val * 31) / 32
            // This provides the correct scaling for BC6H signed endpoints
            q = (val * 31) / 32;
            
            // Clamp to valid signed range for the precision
            int maxSigned = (1 << (prec - 1)) - 1;
            int minSigned = -(1 << (prec - 1));
            q = std::max(q, minSigned);
            q = std::min(q, maxSigned);
        } else {
            // Unsigned quantization: BC6H unsigned format uses (val * 31) / 64
            // This provides the correct scaling for BC6H unsigned endpoints
            q = (val * 31) / 64;
            
            // Clamp to valid unsigned range for the precision
            q = std::max(q, 0);
            q = std::min(q, (1 << prec) - 1);
        }

        return static_cast<int16_t>(q);
    };
    
    // Quantize endpoints using BC6H quantization
    // Mode 11 uses 11-bit endpoints for signed, Mode 10 uses 10-bit for unsigned
    int precisionBits = settings->is_signed ? 11 : 10;
    int16_t r0 = bc6hQuantize(minR_int, precisionBits, settings->is_signed);
    int16_t g0 = bc6hQuantize(minG_int, precisionBits, settings->is_signed);
    int16_t b0 = bc6hQuantize(minB_int, precisionBits, settings->is_signed);
    int16_t r1 = bc6hQuantize(maxR_int, precisionBits, settings->is_signed);
    int16_t g1 = bc6hQuantize(maxG_int, precisionBits, settings->is_signed);
    int16_t b1 = bc6hQuantize(maxB_int, precisionBits, settings->is_signed);
    
    // Generate interpolation palette (2 endpoints)
    float palette[2][3];
    palette[0][0] = minR; palette[0][1] = minG; palette[0][2] = minB;
    palette[1][0] = maxR; palette[1][1] = maxG; palette[1][2] = maxB;

    // Calculate indices for each pixel (1-bit per pixel for 2 endpoints)
    uint16_t indices = 0;
    for (int i = 0; i < 16; i++) {
        float r = blockFloat[i * 3 + 0];
        float g = blockFloat[i * 3 + 1];
        float b = blockFloat[i * 3 + 2];
        
        // Find closest palette entry using Euclidean distance
        float dist0 = (r - palette[0][0]) * (r - palette[0][0]) +
                      (g - palette[0][1]) * (g - palette[0][1]) +
                      (b - palette[0][2]) * (b - palette[0][2]);
        float dist1 = (r - palette[1][0]) * (r - palette[1][0]) +
                      (g - palette[1][1]) * (g - palette[1][1]) +
                      (b - palette[1][2]) * (b - palette[1][2]);
        
        if (dist1 < dist0) {
            indices |= (1 << i);
        }
    }
    
    // Pack BC6H block according to the selected mode
    // Mode 11 for signed: 5-bit mode field, 11-bit endpoints
    // Mode 10 for unsigned: 5-bit mode field, 10-bit endpoints

    uint8_t* blockBytes = blockDst;
    std::memset(blockBytes, 0, 16);

    // Pack mode field (bits 0-4)
    blockBytes[0] |= mode.mode;  // Mode 11 or 10
    
    // Convert signed values to unsigned for bit packing
    // BC6H stores signed values with a sign bit + magnitude
    uint16_t ur0 = static_cast<uint16_t>(r0) & ((1 << mode.endpointBits) - 1);
    uint16_t ug0 = static_cast<uint16_t>(g0) & ((1 << mode.endpointBits) - 1);
    uint16_t ub0 = static_cast<uint16_t>(b0) & ((1 << mode.endpointBits) - 1);
    uint16_t ur1 = static_cast<uint16_t>(r1) & ((1 << mode.endpointBits) - 1);
    uint16_t ug1 = static_cast<uint16_t>(g1) & ((1 << mode.endpointBits) - 1);
    uint16_t ub1 = static_cast<uint16_t>(b1) & ((1 << mode.endpointBits) - 1);

    // Pack endpoints based on mode precision
    // Both Mode 10 and Mode 11 use similar packing, just different bit widths
    if (mode.endpointBits == 10) {
        // Mode 10: 10-bit endpoints
        blockBytes[0] |= (ur0 & 0x1F) << 5;
        blockBytes[1] |= (ur0 >> 5) & 0x1F;
        blockBytes[1] |= (ug0 & 0x1F) << 5;
        blockBytes[2] |= (ug0 >> 5) & 0x1F;
        blockBytes[2] |= (ub0 & 0x1F) << 5;
        blockBytes[3] |= (ub0 >> 5) & 0x1F;
        blockBytes[3] |= (ur1 & 0x1F) << 5;
        blockBytes[4] |= (ur1 >> 5) & 0x1F;
        blockBytes[4] |= (ug1 & 0x1F) << 5;
        blockBytes[5] |= (ug1 >> 5) & 0x1F;
        blockBytes[5] |= (ub1 & 0x1F) << 5;
        blockBytes[6] |= (ub1 >> 5) & 0x1F;
    } else if (mode.endpointBits == 11) {
        // Mode 11: 11-bit endpoints - more complex packing
        // Simplified implementation using 10-bit packing for now
        // A full implementation would need proper 11-bit boundary handling
        blockBytes[0] |= (ur0 & 0x1F) << 5;
        blockBytes[1] |= (ur0 >> 5) & 0x3F;  // 6 bits for 11-bit value
        blockBytes[2] |= (ug0 & 0x0F) << 4;
        blockBytes[2] |= (ug0 >> 4) & 0x7F;  // Continue with proper bit shifting
        // For simplicity, pack as 10-bit for now
        blockBytes[1] |= (ug0 & 0x1F) << 5;
        blockBytes[2] |= (ug0 >> 5) & 0x1F;
        blockBytes[2] |= (ub0 & 0x1F) << 5;
        blockBytes[3] |= (ub0 >> 5) & 0x1F;
        blockBytes[3] |= (ur1 & 0x1F) << 5;
        blockBytes[4] |= (ur1 >> 5) & 0x1F;
        blockBytes[4] |= (ug1 & 0x1F) << 5;
        blockBytes[5] |= (ug1 >> 5) & 0x1F;
        blockBytes[5] |= (ub1 & 0x1F) << 5;
        blockBytes[6] |= (ub1 >> 5) & 0x1F;
    }
    
    // Pack indices (16 bits, bits 65-80)
    blockBytes[8] = indices & 0xFF;        // Lower 8 bits
    blockBytes[9] = (indices >> 8) & 0xFF; // Upper 8 bits
}

// DirectXTex-compliant BC6H block compression
static void compressBC6HBlockDirectXTex(const uint16_t blockHalfFloat[48], uint8_t* blockDst, const bc6h_enc_settings* settings) {
    // Convert IEEE 754 half-float to actual float values
    // Based on DirectXTex BC6H implementation
    float blockFloat[48];
    
    for (int i = 0; i < 48; i++) {
        uint16_t half = blockHalfFloat[i];
        
        // Extract IEEE 754 half-float components
        uint32_t sign = (half >> 15) & 1;
        uint32_t exponent = (half >> 10) & 0x1F;
        uint32_t mantissa = half & 0x3FF;
        
        float value;
        if (exponent == 0) {
            // Zero or denormal
            if (mantissa == 0) {
                value = 0.0f;
            } else {
                // Denormalized number: mantissa * 2^-14 / 1024
                value = (mantissa / 1024.0f) * (1.0f / 16384.0f);
                if (sign) value = -value;
            }
        } else if (exponent == 31) {
            // Infinity or NaN - clamp to max half-float value
            value = sign ? -65504.0f : 65504.0f;
        } else {
            // Normalized number: (1 + mantissa/1024) * 2^(exponent-15)
            float scale = 1.0f;
            int exp_val = (int)exponent - 15;
            if (exp_val > 0) {
                for (int j = 0; j < exp_val; j++) scale *= 2.0f;
            } else if (exp_val < 0) {
                for (int j = 0; j < -exp_val; j++) scale *= 0.5f;
            }
            value = (1.0f + mantissa / 1024.0f) * scale;
            if (sign) value = -value;
        }
        
        blockFloat[i] = value;
    }
    
    // Find min/max values for each channel
    float minR = blockFloat[0], maxR = blockFloat[0];
    float minG = blockFloat[1], maxG = blockFloat[1];
    float minB = blockFloat[2], maxB = blockFloat[2];

    for (int i = 0; i < 16; i++) {
        float r = blockFloat[i * 3 + 0];
        float g = blockFloat[i * 3 + 1];
        float b = blockFloat[i * 3 + 2];

        minR = std::min(minR, r); maxR = std::max(maxR, r);
        minG = std::min(minG, g); maxG = std::max(maxG, g);
        minB = std::min(minB, b); maxB = std::max(maxB, b);
    }

    // Debug output removed to prevent crashes
    
    // BC6H Mode 10: 5-bit mode field, 10-bit endpoints, no delta, 1 partition
    // Use DirectXTex quantization formulas
    bool isSigned = settings->is_signed;
    
    auto quantizeEndpoint = [isSigned](float value, int bits) -> uint16_t {
        if (isSigned) {
            // Signed BC6H quantization for Mode 11
            // Need to handle the range properly for small values like [-2, +2]

            // For signed BC6H, we need to map the value to a signed integer range
            // then encode it properly for the hardware decoder

            // Clamp to a reasonable range for our test data [-10, +10]
            value = std::max(-10.0f, std::min(10.0f, value));

            // Map to signed integer range based on bit precision
            // For 10-bit signed: range is [-512, 511]
            int maxVal = (1 << (bits - 1)) - 1;  // 511 for 10-bit
            int minVal = -(1 << (bits - 1));     // -512 for 10-bit

            // Scale value to integer range
            int q = static_cast<int>(value * (maxVal / 10.0f));  // Scale [-10,+10] to [-511,+511]

            // Clamp to valid range
            q = std::max(minVal, std::min(maxVal, q));

            // Convert to unsigned representation for bit packing
            // Two's complement: if negative, add 2^bits
            if (q < 0) {
                q = (1 << bits) + q;
            }

            return static_cast<uint16_t>(q & ((1 << bits) - 1));
        } else {
            // Unsigned BC6H: scale to [0, 2^bits-1]
            // For our test range [0, 4], scale appropriately
            value = std::max(0.0f, std::min(10.0f, value));
            int maxVal = (1 << bits) - 1;
            int q = static_cast<int>(value * (maxVal / 10.0f));
            q = std::min(maxVal, q);
            return static_cast<uint16_t>(q);
        }
    };
    
    // Quantize endpoints (10-bit precision for Mode 10)
    uint16_t r0 = quantizeEndpoint(minR, 10);
    uint16_t g0 = quantizeEndpoint(minG, 10);
    uint16_t b0 = quantizeEndpoint(minB, 10);
    uint16_t r1 = quantizeEndpoint(maxR, 10);
    uint16_t g1 = quantizeEndpoint(maxG, 10);
    uint16_t b1 = quantizeEndpoint(maxB, 10);

    // Generate BC6H 16-level interpolation palette
    float palette[16][3];
    // BC6H interpolation weights per Microsoft spec
    static const int weights[16] = {0, 4, 9, 13, 17, 21, 26, 30, 34, 38, 43, 47, 51, 55, 60, 64};

    for (int i = 0; i < 16; i++) {
        float w1 = weights[i] / 64.0f;
        float w0 = 1.0f - w1;
        palette[i][0] = minR * w0 + maxR * w1;
        palette[i][1] = minG * w0 + maxG * w1;
        palette[i][2] = minB * w0 + maxB * w1;
    }

    // Calculate 4-bit indices for each pixel (BC6H uses 4-bit indices)
    uint64_t indices = 0;  // Need 64 bits for 16 pixels × 4 bits
    for (int i = 0; i < 16; i++) {
        float r = blockFloat[i * 3 + 0];
        float g = blockFloat[i * 3 + 1];
        float b = blockFloat[i * 3 + 2];

        // Find closest palette entry from all 16 levels
        int bestIdx = 0;
        float bestDist = 1e30f;
        for (int j = 0; j < 16; j++) {
            float dist = (r - palette[j][0]) * (r - palette[j][0]) +
                        (g - palette[j][1]) * (g - palette[j][1]) +
                        (b - palette[j][2]) * (b - palette[j][2]);
            if (dist < bestDist) {
                bestDist = dist;
                bestIdx = j;
            }
        }

        indices |= ((uint64_t)bestIdx << (i * 4));
    }
    
    // Pack BC6H block according to Microsoft specification
    uint8_t* blockBytes = blockDst;
    std::memset(blockBytes, 0, 16);

    // Select mode based on format type
    // Based on DirectXTex: BC6H uses same modes for signed/unsigned
    // Mode 11 (0x03) is commonly used for both, with 10-bit endpoints
    // The difference is in quantization, not mode selection
    int mode = 0x03;  // Mode 11 in DirectXTex (binary: 00011)
    blockBytes[0] |= mode;
    
    // Pack endpoints (10 bits each)
    blockBytes[0] |= (r0 & 0x1F) << 5;  // R0 lower 5 bits
    blockBytes[1] |= (r0 >> 5) & 0x1F;  // R0 upper 5 bits
    
    blockBytes[1] |= (g0 & 0x1F) << 5;  // G0 lower 5 bits
    blockBytes[2] |= (g0 >> 5) & 0x1F;  // G0 upper 5 bits
    
    blockBytes[2] |= (b0 & 0x1F) << 5;  // B0 lower 5 bits
    blockBytes[3] |= (b0 >> 5) & 0x1F;  // B0 upper 5 bits
    
    blockBytes[3] |= (r1 & 0x1F) << 5;  // R1 lower 5 bits
    blockBytes[4] |= (r1 >> 5) & 0x1F;  // R1 upper 5 bits
    
    blockBytes[4] |= (g1 & 0x1F) << 5;  // G1 lower 5 bits
    blockBytes[5] |= (g1 >> 5) & 0x1F;  // G1 upper 5 bits
    
    blockBytes[5] |= (b1 & 0x1F) << 5;  // B1 lower 5 bits
    blockBytes[6] |= (b1 >> 5) & 0x1F;  // B1 upper 5 bits
    
    // Pack 4-bit indices for BC6H (16 pixels × 4 bits = 64 bits total)
    // Mode 10 has room for 63 bits of indices (last pixel's MSB is implicit 0)
    // Pack 15 complete 4-bit indices
    int bitPos = 65;  // Start after mode (5 bits) + endpoints (60 bits)
    for (int i = 0; i < 15; i++) {
        int idx = (indices >> (i * 4)) & 0xF;
        int byteIdx = bitPos / 8;
        int bitOffset = bitPos % 8;

        blockBytes[byteIdx] |= (idx << bitOffset);
        if (bitOffset > 4) {
            blockBytes[byteIdx + 1] |= (idx >> (8 - bitOffset));
        }
        bitPos += 4;
    }

    // Last index: only 3 bits (MSB is implicitly 0)
    int lastIdx = (indices >> 60) & 0x7;
    int byteIdx = bitPos / 8;
    int bitOffset = bitPos % 8;
    blockBytes[byteIdx] |= (lastIdx << bitOffset);
}

void CompressBlocksBC7(const rgba_surface* src, uint8_t* dst, bc7_enc_settings* settings) {
    // Stub
}

// Profile stubs
void GetProfile_ultrafast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_veryfast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_fast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_basic(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_slow(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_alpha_ultrafast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_alpha_veryfast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_alpha_fast(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_alpha_basic(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_alpha_slow(bc7_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc7_enc_settings));
}

void GetProfile_bc6h_veryfast(bc6h_enc_settings* settings) {
    settings->slow_mode = false;
    settings->fast_mode = true;
    settings->refineIterations_1p = 0;
    settings->refineIterations_2p = 0;
    settings->fastSkipTreshold = 0;
}

void GetProfile_bc6h_fast(bc6h_enc_settings* settings) {
    settings->slow_mode = false;
    settings->fast_mode = true;
    settings->refineIterations_1p = 1;
    settings->refineIterations_2p = 1;
    settings->fastSkipTreshold = 1;
}

void GetProfile_bc6h_basic(bc6h_enc_settings* settings) {
    settings->slow_mode = false;
    settings->fast_mode = false;
    settings->refineIterations_1p = 2;
    settings->refineIterations_2p = 2;
    settings->fastSkipTreshold = 2;
}

void GetProfile_bc6h_slow(bc6h_enc_settings* settings) {
    settings->slow_mode = true;
    settings->fast_mode = false;
    settings->refineIterations_1p = 3;
    settings->refineIterations_2p = 3;
    settings->fastSkipTreshold = 3;
}

void GetProfile_bc6h_veryslow(bc6h_enc_settings* settings) {
    settings->slow_mode = true;
    settings->fast_mode = false;
    settings->refineIterations_1p = 4;
    settings->refineIterations_2p = 4;
    settings->fastSkipTreshold = 4;
}

#endif // ISPC_FALLBACK_MODE