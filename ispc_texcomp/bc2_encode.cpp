// Minimal BC2 encoder implementation
// BC2 = 64 bits of explicit 4-bit alpha + 64 bits of BC1 color

#include "ispc_texcomp.h"
#include <cstring>
#include <algorithm>

extern "C" void CompressBlocksBC1(const rgba_surface* src, uint8_t* dst);

extern "C" void CompressBlocksBC2(const rgba_surface* src, uint8_t* dst)
{
    // BC2 block layout:
    // Bytes 0-7:   Alpha block (4 bits per pixel, 16 pixels)
    // Bytes 8-15:  Color block (BC1 format)

    const int block_width = 4;
    const int block_height = 4;
    const int blocks_x = (src->width + 3) / 4;
    const int blocks_y = (src->height + 3) / 4;

    for (int block_y = 0; block_y < blocks_y; block_y++) {
        for (int block_x = 0; block_x < blocks_x; block_x++) {
            uint8_t* block_dst = dst + (block_y * blocks_x + block_x) * 16;

            // First 8 bytes: 4-bit alpha values
            // BC2 stores 4 bits per alpha, packed as 2 pixels per byte
            for (int y = 0; y < block_height; y++) {
                uint16_t alpha_row = 0;
                for (int x = 0; x < block_width; x++) {
                    int px = block_x * block_width + x;
                    int py = block_y * block_height + y;

                    // Get alpha value and clamp to image bounds
                    uint8_t alpha = 255;
                    if (px < src->width && py < src->height) {
                        const uint8_t* pixel = src->ptr + py * src->stride + px * 4;
                        alpha = pixel[3];
                    }

                    // Quantize to 4 bits and pack
                    uint8_t alpha_4bit = (alpha >> 4) & 0x0F;
                    alpha_row |= (alpha_4bit << (x * 4));
                }
                // Store as little-endian
                block_dst[y * 2] = alpha_row & 0xFF;
                block_dst[y * 2 + 1] = (alpha_row >> 8) & 0xFF;
            }

            // Last 8 bytes: BC1 color block
            // Create a temporary surface for this block
            rgba_surface block_surface;
            uint8_t block_data[64]; // 4x4 pixels * 4 bytes

            // Copy block data
            for (int y = 0; y < block_height; y++) {
                for (int x = 0; x < block_width; x++) {
                    int px = block_x * block_width + x;
                    int py = block_y * block_height + y;

                    uint8_t* dst_pixel = block_data + (y * block_width + x) * 4;

                    if (px < src->width && py < src->height) {
                        const uint8_t* src_pixel = src->ptr + py * src->stride + px * 4;
                        memcpy(dst_pixel, src_pixel, 4);
                    } else {
                        // Fill with black for out-of-bounds pixels
                        dst_pixel[0] = 0;
                        dst_pixel[1] = 0;
                        dst_pixel[2] = 0;
                        dst_pixel[3] = 255;
                    }
                }
            }

            block_surface.ptr = block_data;
            block_surface.width = block_width;
            block_surface.height = block_height;
            block_surface.stride = block_width * 4;

            // Compress color using BC1
            CompressBlocksBC1(&block_surface, block_dst + 8);
        }
    }
}