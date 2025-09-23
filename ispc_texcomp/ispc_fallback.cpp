// Fallback stubs when ISPC compiler is not available
#ifdef ISPC_FALLBACK_MODE

#include "ispc_texcomp.h"
#include <cstring>
#include <vector>
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

void CompressBlocksBC3(const rgba_surface* src, uint8_t* dst) {
#ifdef ISPC_FALLBACK_MODE
    // Use the simple BC3 encoder from bc_compression.h
    std::vector<uint8_t> compressed = bc::compressToBC(
        VK_FORMAT_BC3_UNORM_BLOCK, src->ptr, src->width, src->height);
    if (!compressed.empty()) {
        memcpy(dst, compressed.data(), compressed.size());
    }
#endif
}

void CompressBlocksBC4(const rgba_surface* src, uint8_t* dst) {
    // Stub
}

void CompressBlocksBC5(const rgba_surface* src, uint8_t* dst) {
    // Stub
}

void CompressBlocksBC6H(const rgba_surface* src, uint8_t* dst, bc6h_enc_settings* settings) {
    // Stub
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
    std::memset(settings, 0, sizeof(bc6h_enc_settings));
}

void GetProfile_bc6h_fast(bc6h_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc6h_enc_settings));
}

void GetProfile_bc6h_basic(bc6h_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc6h_enc_settings));
}

void GetProfile_bc6h_slow(bc6h_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc6h_enc_settings));
}

void GetProfile_bc6h_veryslow(bc6h_enc_settings* settings) {
    std::memset(settings, 0, sizeof(bc6h_enc_settings));
}

#endif // ISPC_FALLBACK_MODE