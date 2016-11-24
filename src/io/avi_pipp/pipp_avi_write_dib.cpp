#include "pipp_avi_write_dib.h"
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <sstream>

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/param.h>		// define or not BSD macro
#endif

// 64-bit fseek for various platforms
#ifdef __linux__
#define fseek64 fseeko64  // Linux
#define ftell64 ftello64  // Linux
#elif defined (__APPLE__) && defined(__MACH__)
#define fseek64 fseeko  // OS X
#define ftell64 ftello  // OS X
#elif defined(BSD)
#define fseek64 fseeko  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#define ftell64 ftello  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#else
#define fseek64 _fseeki64  // Windows
#define ftell64 _ftelli64  // Windows
#endif

#define MONO_DATA_TIGHTLY_PACKED 1


// ------------------------------------------
// Constructor
// ------------------------------------------
c_pipp_avi_write_dib::c_pipp_avi_write_dib() :
  m_line_gap(0)
{
    m_vids_stream_header.handler.u32 = FCC_DIB;  // Override value from base class
    m_bitmap_info_header.compression.u32 = 0;  // No compression
}


// ------------------------------------------
// Set codec specific values
// ------------------------------------------
int32_t c_pipp_avi_write_dib::set_codec_values()
{
    // Set flag to write colour table if required
    if (!m_colour) {
        m_write_colour_table = 1;
    }

    // Calculate line gap
    m_line_gap = (m_bytes_per_pixel * m_width) % 4;
    if (m_line_gap != 0) {
        m_line_gap = 4 - m_line_gap;
    }

    // Calculate the frame size
    m_frame_size = (m_width * m_bytes_per_pixel + m_line_gap) * m_height;

    return 0;
}


// ------------------------------------------
// Write frame to AVI file
// ------------------------------------------
bool c_pipp_avi_write_dib::write_frame(
    uint8_t  *data,
    int32_t colour,
    uint32_t bpp,
    void *extra_data)
{
    // Remove unused argument warnings
    (void)extra_data;

    // Early return if no file is open
    if (!m_open) {
        return true;
    }
    
    if (colour < 0 || colour > 2) {
        colour = 0;
    }

    // Indicate that a frame is about to be added
    frame_added();
    fwrite_error_check(&m_00db_chunk_header, 1, sizeof(m_00db_chunk_header), mp_avi_file);

    uint8_t *buffer;
    int32_t line_length = m_width * m_bytes_per_pixel;
    if (m_line_gap == 0 && m_bytes_per_pixel == 3 && bpp == 1) {
        // Use supplied data directly
        buffer = data;
    } else {
        // Create version of image with line gaps in
        buffer = m_temp_buffer.get_buffer((m_width * m_bytes_per_pixel + m_line_gap)* m_height);
        
        if (bpp == 1) {
            if (m_bytes_per_pixel == 3) {
                // Colour version 
                for (int32_t y = 0; y < m_height; y++) {
                    uint8_t *src_ptr = data + (y * line_length);
                    uint8_t *dst_ptr = buffer + (y * (line_length + m_line_gap));
                    memcpy(dst_ptr, src_ptr, line_length);
                }
            } else {
                // Mono version
                uint8_t *src_ptr = data + colour;
                uint8_t *dst_ptr = buffer;
                for (int32_t y = 0; y < m_height; y++) {
                    for (int32_t x = 0; x < m_width; x++) {
                        *dst_ptr++ = *src_ptr;
#ifdef MONO_DATA_TIGHTLY_PACKED
                        src_ptr++;
#else
                        src_ptr += 3;
#endif
                    }

                    dst_ptr += m_line_gap;
                }
            }
        } else {  // Bytes per sample == 2
            if (m_bytes_per_pixel == 3) {
                // Colour version 
                for (int32_t y = 0; y < m_height; y++) {
                    uint16_t *src_ptr = (uint16_t *)data + (y * line_length);
                    uint8_t *dst_ptr = buffer + (y * (line_length + m_line_gap));
                    for (int32_t x = 0; x < line_length; x++) {
                        *dst_ptr++ = *src_ptr++ >> 8;
                    }
                }
            } else {
                // Mono version
                uint16_t *src_ptr = (uint16_t *)data + colour;
                uint8_t *dst_ptr = buffer;
                for (int32_t y = 0; y < m_height; y++) {
                    for (int32_t x = 0; x < m_width; x++) {
                        *dst_ptr++ = *src_ptr >> 8;
#ifdef MONO_DATA_TIGHTLY_PACKED
                        src_ptr++;
#else
                        src_ptr += 3;
#endif
                    }

                    dst_ptr += m_line_gap;
                }
            }
        }
    }

    // Write image data to file
    m_last_frame_pos = ftell64(mp_avi_file);  // Grab position of last file
    fwrite_error_check(buffer , 1 , (m_width * m_bytes_per_pixel + m_line_gap) * m_height, mp_avi_file);

    // Tidy up after write failures
    if (m_file_write_error) {
        fclose(mp_avi_file);
        m_open = false;
    }

    bool ret = m_file_write_error;
    m_file_write_error = false;
    return ret;
}
