#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include "pipp_avi_write.h"

#include <cwchar>
#include <memory>
#include <sstream>

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <sys/param.h>		// define or not BSD macro
#endif

// 64-bit fseek for various platforms
#ifdef __linux__
#define fseek64 fseeko64  // Linux
#define ftell64 ftello64  // Linux
#elif defined (__APPLE__) && defined (__MACH__)
#define fseek64 fseeko  // OS X
#define ftell64 ftello  // OS X
#elif defined(BSD)
#define fseek64 fseeko  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#define ftell64 ftello  // DragonFly BSD, FreeBSD, OpenBSD, NetBSD
#elif defined (__FreeBSD_kernel__) && defined (__GLIBC__)
#define fseek64 fseeko  // KFreeBSD
#define ftell64 ftello  // KFreeBSD
#else
#define fseek64 _fseeki64  // Windows
#define ftell64 _ftelli64  // Windows
#endif

using namespace std;


//
// index_type codes
//
#define AVI_INDEX_OF_INDEXES   0x00  // when each entry in aIndex array points to an index chunk
#define AVI_INDEX_OF_CHUNKS    0x01  // when each entry in aIndex array points to a chunk in the file
#define AVI_INDEX_IS_DATA      0x80  // when each entry is aIndex is really the data

//
// index_sub_type codes for INDEX_OF_CHUNKS
//
#define AVI_INDEX_2FIELD 0x01 // when fields within frames are also indexed


// ------------------------------------------
// Constructor
// ------------------------------------------
c_pipp_avi_write::c_pipp_avi_write() :
    mp_avi_file(NULL),
    m_open(false),
    m_split_count(0),
    m_old_avi_format(0),
    m_width(0),
    m_height(0),
    m_frame_size(0),
    m_colour(false),
    m_write_colour_table(false),
    m_total_frame_count(0),
    m_current_frame_count(0),
    m_riff_count(0),
    m_bytes_per_pixel(1),
    m_last_frame_pos(0),
    m_file_write_error(false)
{
    // Initialise RIFF AVI header
    m_avi_riff_header.list.u32 = FCC_RIFF;
    m_avi_riff_header.size = 0;
    m_avi_riff_header.four_cc.u32 = FCC_AVI;

    // Initialise the RIFF AVIX header
    m_avix_riff_header.list.u32 = FCC_RIFF;
    m_avix_riff_header.size = 0;
    m_avix_riff_header.four_cc.u32 = FCC_AVIX;

    // Initialise hdrl list header
    m_hdrl_list_header.list.u32 = FCC_LIST;
    m_hdrl_list_header.size = 0;
    m_hdrl_list_header.four_cc.u32 = FCC_hdrl;

    // Initialise avih chunk header
    m_avih_chunk_header.four_cc.u32 = FCC_avih;
    m_avih_chunk_header.size = sizeof(s_main_avi_header);

    // Initialise the junk chunk header
    m_junk_chunk_header.four_cc.u32 = FCC_JUNK;
    m_junk_chunk_header.size = 0;

    // Initialise the main AVI header
    m_main_avih_header.micro_sec_per_frame = 100000;  // 10 fps by default
    m_main_avih_header.max_bytes_per_sec = 0;
    m_main_avih_header.padding_granularity = 0;
    m_main_avih_header.flags = AVIF_HASINDEX;// | AVIF_MUSTUSEINDEX;
    m_main_avih_header.total_frames = 0;  // Increment as frames are added
    m_main_avih_header.initial_frames = 0;  // Always 0
    m_main_avih_header.streams = 1;  // Always 1 stream
    m_main_avih_header.suggested_buffer_size = 0;  // Fill in later
    m_main_avih_header.width = 0;  // Fill in later
    m_main_avih_header.height = 0;  // Fill in later
    m_main_avih_header.reserved[0] = 0;  // Always 0
    m_main_avih_header.reserved[1] = 0;  // Always 0
    m_main_avih_header.reserved[2] = 0;  // Always 0
    m_main_avih_header.reserved[3] = 0;  // Always 0

    // Initialise strl list header
    m_strl_list_header.list.u32 = FCC_LIST;
    m_strl_list_header.size = 0;
    m_strl_list_header.four_cc.u32 = FCC_strl;

    // Initialise strh chunk header
    m_strh_chunk_header.four_cc.u32 = FCC_strh;
    m_strh_chunk_header.size = sizeof(s_avi_stream_header);

    // Initialise vids stream header
    m_vids_stream_header.type.u32 = FCC_vids;
    m_vids_stream_header.handler.u32 = FCC_DIB;
    m_vids_stream_header.flags = 0;
    m_vids_stream_header.priority = 0;
    m_vids_stream_header.language = 0;
    m_vids_stream_header.initial_frames = 0;
    m_vids_stream_header.scale = 1;
    m_vids_stream_header.rate = 10;  // 10 fps by default
    m_vids_stream_header.start = 0;
    m_vids_stream_header.length = 0;  // Fill in later
    m_vids_stream_header.suggested_buffer_size = 0;  // Fill in later
    m_vids_stream_header.quality = 0xFFFFFFFF;
    m_vids_stream_header.sample_size = 0;
    m_vids_stream_header.frame.left = 0;
    m_vids_stream_header.frame.top = 0;
    m_vids_stream_header.frame.right = 0;
    m_vids_stream_header.frame.bottom = 0;

    // Initialise STRF chunk header
    m_strf_chunk_header.four_cc.u32 = FCC_strf;
    m_strf_chunk_header.size = 0;  // Fill in later

    // Initialise the BITMAP info header
    m_bitmap_info_header.size = 0;  // Fill in later
    m_bitmap_info_header.width = 0;  // Fill in later
    m_bitmap_info_header.height = 0;  // Fill in later
    m_bitmap_info_header.planes = 1;
    m_bitmap_info_header.bit_count = 8;  // Fill in later
    m_bitmap_info_header.compression.u32 = 0;  // No compression
    m_bitmap_info_header.size_image = 0;  // Fill in later
    m_bitmap_info_header.x_pels_per_meter = 0;
    m_bitmap_info_header.y_pels_per_meter = 0;
    m_bitmap_info_header.clr_used = 0;  // Fill in later
    m_bitmap_info_header.clr_important = 0;

    // Initialise the indx chunk header
    m_indx_chunk_header.four_cc.u32 = FCC_indx;
    m_indx_chunk_header.size = sizeof(s_avi_superindex_header)
                           + sizeof(s_avi_superindex_entry) * NUMBER_SUPERINDEX_ENTRIES;

    // Initialise the AVI superindex header
    m_avi_superindex_header.longs_per_entry = 4;
    m_avi_superindex_header.index_sub_type = 0;//AVI_INDEX_OF_CHUNKS;
    m_avi_superindex_header.index_type = AVI_INDEX_OF_INDEXES;
    m_avi_superindex_header.entries_in_use = 0;  // Will be increment as needed
    m_avi_superindex_header.chunk_id.u32 = FCC_00db;
    m_avi_superindex_header.reserved[0] = 0;
    m_avi_superindex_header.reserved[1] = 0;
    m_avi_superindex_header.reserved[2] = 0;

    // Initialise the AVI superindex entries
    for (int x = 0; x < NUMBER_SUPERINDEX_ENTRIES; x++) {
        m_avi_superindex_entries[x].offset = 0;    // Fill in later
        m_avi_superindex_entries[x].size = 0;      // Fill in later - the size of the standard index this entry point to
        m_avi_superindex_entries[x].duration = 0;  // Fill in later - the number of frames in the sub-index
    }

    // Initialise the odml list header
    m_odml_list_header.list.u32 = FCC_LIST;
    m_odml_list_header.size = sizeof(m_odml_list_header.four_cc.u32) + sizeof(m_dmlh_chunk_header) + sizeof(m_extended_avi_header);
    m_odml_list_header.four_cc.u32 = FCC_odml;

    // Initialise the dmlh chunk header
    m_dmlh_chunk_header.four_cc.u32 = FCC_dmlh;
    m_dmlh_chunk_header.size = sizeof(m_extended_avi_header);

    // Initialise the extended AVI header
    m_extended_avi_header.total_frames = 0;  // Fill in as the file closes

    // Initialise the movie list header (AVI RIFF version)
    m_movi_list_header.list.u32 = FCC_LIST;
    m_movi_list_header.size = 0;  // Fill in later
    m_movi_list_header.four_cc.u32 = FCC_movi;

    // Initialise the movie list header (AVIX RIFF version)
    m_movi_avix_list_header.list.u32 = FCC_LIST;
    m_movi_avix_list_header.size = 0;  // Fill in later
    m_movi_avix_list_header.four_cc.u32 = FCC_movi;

    // Initialise 00db chunk header
    m_00db_chunk_header.four_cc.u32 = FCC_00db;
    m_00db_chunk_header.size = 0;  // Fill in later

    // Initialise idx1 chunk header
    m_idx1_chunk_header.four_cc.u32 = FCC_idx1;
    m_idx1_chunk_header.size = 0;  // Fill in with frames * sizeof(s_avi_old_index_entry)

    // Initialise the ix00 chunk header
    m_ix00_chunk_header.four_cc.u32 = FCC_ix00;
    m_ix00_chunk_header.size = 0;  // Fill in later

    // Initialise the standard index header
    m_avi_stdindex_header.longs_per_entry = 2;
    m_avi_stdindex_header.index_sub_type = 0;
    m_avi_stdindex_header.index_type = AVI_INDEX_OF_CHUNKS;
    m_avi_stdindex_header.entries_in_use = 0;
    m_avi_stdindex_header.chunk_id.u32 = FCC_00db;
    m_avi_stdindex_header.base_offset[0] = 0;
    m_avi_stdindex_header.base_offset[1] = 0;
    m_avi_stdindex_header.reserved3 = 0;

    // Initialise the standard index entry
    m_avi_stdindex_entry.offset = 0;
    m_avi_stdindex_entry.size = 0;

    // Initialise avi index entry
    m_avi_index_entry.chunk_id.u32 = FCC_00db;
    m_avi_index_entry.flags = AVIIF_KEYFRAME;
    m_avi_index_entry.offset = 0;  // Fill in as required
    m_avi_index_entry.size = 0;  // Fill in later
}


// ------------------------------------------
// fwrite() function with error checking
// ------------------------------------------
void c_pipp_avi_write::fwrite_error_check(
    const void *ptr,
    size_t size,
    size_t count,
    FILE *p_stream)
{
    if (!m_file_write_error) {  // Do not continue writing after an error has occured
        size_t size_written = fwrite(ptr, size, count, p_stream);
        if (size_written != count) {
            m_file_write_error = true;
        }
    }
}


// ------------------------------------------
// Write headers to file
// ------------------------------------------
void c_pipp_avi_write::write_headers()
{
    if (m_file_write_error) {
        // Early exit if there has already been a file write error
        return;
    }

    // Write RIF Header to file
    fwrite_error_check(&m_avi_riff_header , 1 , sizeof(m_avi_riff_header) , mp_avi_file);

    // Write hdrl list header to file
    fwrite_error_check(&m_hdrl_list_header , 1 , sizeof(m_hdrl_list_header) , mp_avi_file);

    // Write avih chunk header to file
    fwrite_error_check(&m_avih_chunk_header , 1 , sizeof(m_avih_chunk_header) , mp_avi_file);

    // Write Main avi header to file
    fwrite_error_check(&m_main_avih_header , 1 , sizeof(m_main_avih_header) , mp_avi_file);

    // Write strl list header to file
    fwrite_error_check(&m_strl_list_header , 1 , sizeof(m_strl_list_header) , mp_avi_file);

    // Write strh chunk header to file
    fwrite_error_check(&m_strh_chunk_header , 1 , sizeof(m_strh_chunk_header) , mp_avi_file);

    // Write video stream header to file
    fwrite_error_check(&m_vids_stream_header , 1 , sizeof(m_vids_stream_header) , mp_avi_file);

    // Write strf chunk header to file
    fwrite_error_check(&m_strf_chunk_header , 1 , sizeof(m_strf_chunk_header) , mp_avi_file);

    // Write BITMAPINFO header to file
    fwrite_error_check(&m_bitmap_info_header , 1 , sizeof(m_bitmap_info_header) , mp_avi_file);

    // Write colour table to file if required (for DIB mono images)
    if (m_write_colour_table) {
        uint8_t colour_table[256 * 4];
        for (int32_t x = 0; x < 256 * 4; x++) {
            if (x % 4 != 3) {
                colour_table[x] = x / 4;
            } else {
                colour_table[x] = 0;
            }
        }

        // Write colour table to file
        fwrite_error_check(colour_table , 1 , 256 * 4, mp_avi_file);
    }

    if (m_old_avi_format != 0) {
        // Junk chunk moves next pos to 0x2000
        m_junk_chunk_header.size = 0x2000 - (int32_t)ftell64(mp_avi_file) - sizeof(m_junk_chunk_header);

        // Write junk header to file
        fwrite_error_check(&m_junk_chunk_header , 1 , sizeof(m_junk_chunk_header) , mp_avi_file);

        // Write junk data to file
        uint8_t junk_byte = 0;
        for (uint32_t junk_count = 0; junk_count < m_junk_chunk_header.size; junk_count++) {
            fwrite_error_check(&junk_byte , 1 , 1, mp_avi_file);
        }
    } else {
        // These fields are not present with the old AVI format

        // Write indx chunk header to file
        fwrite_error_check(&m_indx_chunk_header , 1 , sizeof(m_indx_chunk_header) , mp_avi_file);

        // Write AVI Superindex header to file
        fwrite_error_check(&m_avi_superindex_header , 1 , sizeof(m_avi_superindex_header), mp_avi_file);

        // Write AVI Superindex entries to file
        fwrite_error_check(m_avi_superindex_entries , 1 , sizeof(m_avi_superindex_entries[0]) * NUMBER_SUPERINDEX_ENTRIES, mp_avi_file);

        // Write odml list header to file
        fwrite_error_check(&m_odml_list_header , 1 , sizeof(m_odml_list_header) , mp_avi_file);

        // Write dmlh chunk header to file
        fwrite_error_check(&m_dmlh_chunk_header , 1 , sizeof(m_dmlh_chunk_header) , mp_avi_file);

        // Write extended AVI header to file
        fwrite_error_check(&m_extended_avi_header , 1 , sizeof(m_extended_avi_header) , mp_avi_file);
    }

    // Write movi list header
    fwrite_error_check(&m_movi_list_header , 1 , sizeof(m_movi_list_header) , mp_avi_file);
}


// ------------------------------------------
// A new frame has been added
// ------------------------------------------
void c_pipp_avi_write::frame_added()
{
    // Split AVI files to prevent max size from being exceeded
    if (m_old_avi_format != 0) {
        if (m_current_frame_count == m_max_frames_in_first_riff) {
            split_close();
            split_create();
        }
    } else {
        // Not old format
        if (m_riff_count == (NUMBER_SUPERINDEX_ENTRIES-1) && m_current_frame_count == m_max_frames_in_other_riffs) {
            split_close();
            split_create();
        }
    }

    if (m_total_frame_count == 0) {
        // Grab position of first frame in this RIFF for the base offset
        m_avi_superindex_header.entries_in_use++;
        uint64_t base_offset = ftell64(mp_avi_file) + sizeof(m_00db_chunk_header);
        //uint64_t base_offset = ftell64(avi_fp) - frame_size;
        m_avi_stdindex_header.base_offset[1] = (uint32_t)(base_offset >> 32);
        m_avi_stdindex_header.base_offset[0] = (uint32_t)(base_offset & 0xFFFFFFFF);
    }
    
    m_total_frame_count++;  // Increment frame counts

    if (m_old_avi_format == 0) {
        if ((m_riff_count == 0 && m_current_frame_count == m_max_frames_in_first_riff)
         || (m_riff_count > 0 && m_current_frame_count == m_max_frames_in_other_riffs)) {
            // This frame needs to be in a new RIFF
            finish_riff();  // Finish this RIFF

            // Start the next RIFF
            fwrite_error_check(&m_avix_riff_header , 1, sizeof(m_avix_riff_header), mp_avi_file);

            // Start the new movi LIST
            fwrite_error_check(&m_movi_avix_list_header , 1, sizeof(m_movi_avix_list_header), mp_avi_file);

            // Grab position of first frame in this RIFF for the base offset
            m_avi_superindex_header.entries_in_use++;
            uint64_t base_offset = ftell64(mp_avi_file) + sizeof(m_00db_chunk_header);
            m_avi_stdindex_header.base_offset[1] = (uint32_t)(base_offset >> 32);
            m_avi_stdindex_header.base_offset[0] = (uint32_t)(base_offset & 0xFFFFFFFF);
            m_current_frame_count = 0;
        }
    }

    m_current_frame_count++;
}


// ------------------------------------------
// Create a new AVI file
// ------------------------------------------
bool c_pipp_avi_write::create(
    const char *filename,
    int32_t width,
    int32_t height,
    bool colour,
    int32_t fps_rate,
    int32_t fps_scale,
    int32_t old_avi_format,
    int32_t quality,
    void *extra_data)
{
    // Remove unused argument warnings
    (void)quality;
    (void)extra_data;

    // Set member variables
    m_width = width;
    m_height = height;
    m_colour = colour;
    m_old_avi_format = old_avi_format;
    if (!colour) {
        m_bytes_per_pixel = 1;
    } else {
        m_bytes_per_pixel = 3;
    }

    // Call derived class method to set codec specific values
    set_codec_values();

    // Calculate how many frames fit into each RIFF
    if (old_avi_format != 0) {
        // Calculate how many frames can go into the RIFF
        if (old_avi_format == 4) {
            // Max RIFF size = 4GB
            m_max_frames_in_first_riff = 0xFFFFFFFF;  // Maximum RIFF size (4GB - 1)
        } else {
            // Max RIFF size = 2GB
            m_max_frames_in_first_riff = 0x7FFFFFFF;  // Maximum RIFF size (2GB - 1)
        }
        //max_frames_in_first_riff -= sizeof(avi_riff_header) ;
        //max_frames_in_first_riff -= (sizeof(hdrl_list_header) - sizeof(hdrl_list_header.four_cc) + hdrl_list_header.size);
        //max_frames_in_first_riff -= (sizeof(junk_chunk_header) + junk_chunk_header.size);
        m_max_frames_in_first_riff -= 0x2000;  // Junk field always takes us to 0x2000
        m_max_frames_in_first_riff -= sizeof(m_movi_list_header);
        m_max_frames_in_first_riff -= sizeof(m_idx1_chunk_header);
        m_max_frames_in_first_riff /= (sizeof(m_00db_chunk_header) + m_frame_size + sizeof(m_avi_index_entry));
    } else {
        // Calculate how many frames can go into the first RIFF
        m_max_frames_in_first_riff = 0x3FFFFFFF;  // Maximum RIFF size (1GB - 1)
        m_max_frames_in_first_riff -= sizeof(m_avi_riff_header) ;
        m_max_frames_in_first_riff -= (sizeof(m_hdrl_list_header) - sizeof(m_hdrl_list_header.four_cc) + m_hdrl_list_header.size);
        m_max_frames_in_first_riff -= sizeof(m_movi_list_header);
        m_max_frames_in_first_riff -= (sizeof(m_ix00_chunk_header) + sizeof(m_avi_stdindex_header));
        m_max_frames_in_first_riff -= sizeof(m_idx1_chunk_header);
        m_max_frames_in_first_riff /= (sizeof(m_00db_chunk_header) + m_frame_size + sizeof(m_avi_stdindex_entry) + sizeof(m_avi_index_entry));

        // Calculate how many frames can go into the subsequent RIFFs
        m_max_frames_in_other_riffs = 0x7FFFFFFF;  // Maximum RIFF size (2GB - 1)
        m_max_frames_in_other_riffs -= (sizeof(m_avix_riff_header) + sizeof(m_movi_avix_list_header) + sizeof(m_ix00_chunk_header) + sizeof(m_avi_stdindex_header));
        m_max_frames_in_other_riffs /= (sizeof(m_00db_chunk_header) + m_frame_size + sizeof(m_avi_stdindex_entry));
    }

    m_split_count = 0;

    mp_filename.reset(new char[strlen(filename) + 1]);
    // Copy filename into buffer
    strcpy(mp_filename.get(), filename);

    // Get extension
    char *extension = strrchr(mp_filename.get(), '.');
    if (extension == NULL) {
        // No extension found - create one
        mp_extension.reset(new char[strlen(".avi") + 1]);
        strcpy(mp_extension.get(), ".avi");  // Copy extension
    } else {
        mp_extension.reset(new char[strlen(extension) + 1]);
        strcpy(mp_extension.get(), extension);  // Copy extension
        *extension = 0;  // Remove extension from filename
    }

    // Reset counts
    m_total_frame_count = 0;
    m_current_frame_count = 0;
    m_riff_count = 0;

    // Set flag to write colour table if required
    if (!colour) {
        m_write_colour_table = 1;
    }

    // Update AVI structures
    m_main_avih_header.width = width;
    m_main_avih_header.height = height;
    uint64_t us_per_frame = (uint64_t)1000000 * (uint64_t)fps_scale;
    us_per_frame /= fps_rate;
    m_main_avih_header.micro_sec_per_frame = (uint32_t)us_per_frame;
    m_vids_stream_header.rate = fps_rate;
    m_vids_stream_header.scale = fps_scale;
    m_vids_stream_header.frame.right = width;
    m_vids_stream_header.frame.bottom = height;

    // Reset fields to count frames and indexes
    m_avi_superindex_header.entries_in_use = 0;  // Will be increment as needed
    m_extended_avi_header.total_frames = 0;  // Increment as frames are added

    // Set size of strf chunk
    m_strf_chunk_header.size = sizeof(m_bitmap_info_header);
    if (m_write_colour_table == 1) {
        m_strf_chunk_header.size += 256 * 4;
    }

    // Set size of strl LIST
    m_strl_list_header.size = sizeof(m_strl_list_header.four_cc)
                          + sizeof(m_strh_chunk_header) + m_strh_chunk_header.size
                          + sizeof(m_strf_chunk_header) + m_strf_chunk_header.size;

    if (old_avi_format == 0) {
        m_strl_list_header.size = m_strl_list_header.size
                               + sizeof(m_indx_chunk_header) + m_indx_chunk_header.size
                               + sizeof(m_odml_list_header) - sizeof(m_odml_list_header.four_cc) + m_odml_list_header.size;
    }
        
    // Set the size of hdrl LIST
    m_hdrl_list_header.size = sizeof(m_hdrl_list_header.four_cc)
                          + sizeof(m_avih_chunk_header) + m_avih_chunk_header.size
                          + sizeof(m_strl_list_header) - sizeof(m_strl_list_header.four_cc) + m_strl_list_header.size;

    m_main_avih_header.suggested_buffer_size   = m_frame_size + 8;
    m_vids_stream_header.suggested_buffer_size = m_frame_size + 8;

    m_bitmap_info_header.size = sizeof(m_bitmap_info_header);
    m_bitmap_info_header.width = width;
    m_bitmap_info_header.height = height;
    m_bitmap_info_header.bit_count = 8 * m_bytes_per_pixel;

    // Set colour table length
    if (m_bytes_per_pixel == 1) {
          m_bitmap_info_header.clr_used = 256;
    }

    m_00db_chunk_header.size = m_frame_size;

    // Set up MOVI LIST size for AVIX RIFFs as maximum size
    m_movi_avix_list_header.size = sizeof(m_movi_avix_list_header.four_cc)
                               + m_max_frames_in_other_riffs * sizeof(s_chunk_header)       // 00db chunks
                               + m_max_frames_in_other_riffs * m_frame_size                   // frame data
                               + sizeof(m_ix00_chunk_header)                                // ix00 chunk header  
                               + sizeof(m_avi_stdindex_header)                              // Standard index header
                               + m_max_frames_in_other_riffs * sizeof(m_avi_stdindex_entry);  // Standard index entries

    // Set up AVIX RIFF size as maximum size
    m_avix_riff_header.size = sizeof(m_avix_riff_header.four_cc)
                          + sizeof(m_movi_avix_list_header) - sizeof(m_movi_avix_list_header.four_cc)
                          + m_movi_avix_list_header.size;

    mp_avi_file = fopen(filename, "wb+");

    // Check file opened
    // Return if file did not open
    if (mp_avi_file) {
        m_open = true;
    } else {
        m_open = false;
        m_file_write_error = true;
    }

    // Write headers to file
    write_headers();

    // Handle case where file opened but subsequent write failed
    if (m_open && m_file_write_error) {
        fclose(mp_avi_file);
        m_open = false;
    }

    bool ret = m_file_write_error;
    m_file_write_error = false;
    return ret;
}


// ------------------------------------------
// Create a new AVI file
// ------------------------------------------
void c_pipp_avi_write::split_create()
{
    // Increment split count
    m_split_count++;

    // Reset counts
    m_total_frame_count = 0;
    m_current_frame_count = 0;
    m_riff_count = 0;

    // Reset fields to count frames and indexes
    m_avi_superindex_header.entries_in_use = 0;  // Will be increment as needed
    m_extended_avi_header.total_frames = 0;  // Increment as frames are added

    // Call derived class method to set codec specific values
    set_codec_values();

    // Create split filename
    std::unique_ptr<char[]> p_split_filename(new char[strlen(mp_filename.get()) + 3 + strlen(mp_extension.get()) + 1]);
    if (m_split_count == 0) {
        // No split count to be inserted at end of filename
        sprintf(p_split_filename.get(), "%s%s", mp_filename.get(), mp_extension.get());
    } else {
        // Split count required at end of filename
        sprintf(p_split_filename.get(), "%s_%02d%s", mp_filename.get(), m_split_count, mp_extension.get());
    }

    // Open new file
    mp_avi_file = fopen(p_split_filename.get(), "wb+");

    // Check file opened
    // Return if file did not open
    if (mp_avi_file) {
        m_open = true;
    } else {
        m_file_write_error = true;
        m_open = false;
    }

    // Write headers to file
    write_headers();
}


// ------------------------------------------
// Finish the current RIFF
// ------------------------------------------
void c_pipp_avi_write::finish_riff()
{
    if (m_file_write_error) {
        // Early return for previous file write errors
        return;
    }

    // Add odml indexes
    if (m_old_avi_format == 0) {
        m_ix00_chunk_header.size = sizeof(m_avi_stdindex_header)
                               + m_current_frame_count * sizeof(m_avi_stdindex_entry);

        // Grab position of the ix00 chunk
        m_avi_superindex_entries[m_riff_count].duration = m_current_frame_count;
        m_avi_superindex_entries[m_riff_count].offset = ftell64(mp_avi_file);
        m_avi_superindex_entries[m_riff_count].size = sizeof(m_ix00_chunk_header) + m_ix00_chunk_header.size;

        // Write ix00 chunk header to file
        fwrite_error_check(&m_ix00_chunk_header, 1, sizeof(m_ix00_chunk_header), mp_avi_file);

        // Write AVI standard header to file
        m_avi_stdindex_header.entries_in_use = m_current_frame_count;
        fwrite_error_check(&m_avi_stdindex_header, 1, sizeof(m_avi_stdindex_header), mp_avi_file);

        // Write AVI standard indexes to file
        m_avi_stdindex_entry.size = m_frame_size;
        for (int x = 0; x < m_current_frame_count; x++) {
            // Update entry
            m_avi_stdindex_entry.offset = x * (m_frame_size + sizeof(m_00db_chunk_header));
            fwrite_error_check(&m_avi_stdindex_entry, 1, sizeof(m_avi_stdindex_entry), mp_avi_file);
        }
    }

    // Update final fields
    m_extended_avi_header.total_frames = m_total_frame_count;
    m_vids_stream_header.length = m_total_frame_count;

    if (m_riff_count == 0) {
        // This is the first RIFF - update fields that rely on the first first RIFF

        m_main_avih_header.total_frames = m_current_frame_count;

        m_bitmap_info_header.size_image = m_frame_size;
        m_avi_index_entry.size = m_frame_size;
        m_idx1_chunk_header.size = m_current_frame_count * sizeof(s_avi_old_index_entry);
        m_movi_list_header.size = sizeof(m_movi_list_header.four_cc)
                              + m_current_frame_count * sizeof(s_chunk_header)  // 00db chunks
                              + m_current_frame_count * m_frame_size;             // frame data

        if (m_old_avi_format == 0) {
            // Add ix00 index list size on
            m_movi_list_header.size = m_movi_list_header.size
                                  + sizeof(m_ix00_chunk_header)
                                  + m_ix00_chunk_header.size;
        }

        // Write index chunk header to file
        fwrite_error_check(&m_idx1_chunk_header , 1 , sizeof(m_idx1_chunk_header), mp_avi_file);

        // Write AVI 1.0 index entries to file
        m_avi_index_entry.offset = 0x4;

        // Write all entries
        for (int32_t x = 0; x < m_current_frame_count; x++) {
            fwrite_error_check(&m_avi_index_entry , 1 , sizeof(m_avi_index_entry), mp_avi_file);  // Write entry to file
            m_avi_index_entry.offset += (sizeof(s_chunk_header) + m_frame_size);  // Increment offset
        }

        // Get the filesize
        int64_t filesize = ftell64(mp_avi_file);

        // Update final headers
        m_avi_riff_header.size = (uint32_t)filesize - 8;
    } 

    // Grab start position of the next RIFF
    int64_t riff_end_position = ftell64(mp_avi_file);

    // Processing for subsequent RIFFs
    if (m_riff_count > 0 && m_current_frame_count != m_max_frames_in_other_riffs) {
        // This RIFF must be the last RIFF as it does not have the maximum number of frames in it
        // We need to correct the RIFF and LIST sizes as it is not completely full

        // Go back to the start of this RIFF
        fseek64(mp_avi_file, m_riff_start_position, SEEK_SET);

        // Write RIFF header again with correct length now that we know it
        m_avix_riff_header.size = (int32_t)(riff_end_position - m_riff_start_position) - sizeof(m_avix_riff_header) + sizeof(m_avix_riff_header.four_cc);
        fwrite_error_check(&m_avix_riff_header , 1, sizeof(m_avix_riff_header), mp_avi_file);

        // Write the movi LIST header again with the correct length now that we know it
        m_movi_avix_list_header.size = m_avix_riff_header.size - sizeof(m_movi_avix_list_header);
        fwrite_error_check(&m_movi_avix_list_header , 1 , sizeof(m_movi_avix_list_header) , mp_avi_file);

        // Go back to the end of this RIFF
        fseek64(mp_avi_file, riff_end_position, SEEK_SET);
    }

    // Grab start position of the next RIFF
    m_riff_start_position = riff_end_position;

    // Reset current frame count as this RIFF has been closed
    m_current_frame_count = 0;

    // Increment RIFF count
    m_riff_count++;
}



// ------------------------------------------
// Write header and close AVI file
// ------------------------------------------
bool c_pipp_avi_write::close()
{
    if (m_open) {
        // Finish off this RIFF
        finish_riff();

        // Go back to start of file
        fseek64(mp_avi_file, 0, SEEK_SET);

        // Write the updated headers to the file
        write_headers();

        // Note that the AVI file is closed
        m_open = false;

        fclose(mp_avi_file);
        mp_avi_file = NULL;
    }

    bool ret = m_file_write_error;
    m_file_write_error = false;
    return ret;
}


// ------------------------------------------
// Write header and close AVI file
// ------------------------------------------
void c_pipp_avi_write::split_close()
{
    // Finish off this RIFF
    finish_riff();

    // Go back to start of file
    fseek64(mp_avi_file, 0, SEEK_SET);

    // Write the updated headers to the file
    write_headers();

    fclose(mp_avi_file);
    mp_avi_file = NULL;
}


int32_t c_pipp_avi_write::debug_headers() {
    // Debug output
    printf("frame_count: %d\n\n", m_total_frame_count);

    printf("avi_riff_header.list: %c%c%c%c\n",
           m_avi_riff_header.list.chr[0],
           m_avi_riff_header.list.chr[1],
           m_avi_riff_header.list.chr[2],
           m_avi_riff_header.list.chr[3]);

    printf("avi_riff_header.size: 0x%x\n", m_avi_riff_header.size);

    printf("avi_riff_header.four_cc: %c%c%c%c\n\n",
           m_avi_riff_header.four_cc.chr[0],
           m_avi_riff_header.four_cc.chr[1],
           m_avi_riff_header.four_cc.chr[2],
           m_avi_riff_header.four_cc.chr[3]);

    printf("hdrl_list_header.list: %c%c%c%c\n",
           m_hdrl_list_header.list.chr[0],
           m_hdrl_list_header.list.chr[1],
           m_hdrl_list_header.list.chr[2],
           m_hdrl_list_header.list.chr[3]);

    printf("hdrl_list_header.size: 0x%x\n", m_hdrl_list_header.size);

    printf("hdrl_list_header.four_cc: %c%c%c%c\n\n",
           m_hdrl_list_header.four_cc.chr[0],
           m_hdrl_list_header.four_cc.chr[1],
           m_hdrl_list_header.four_cc.chr[2],
           m_hdrl_list_header.four_cc.chr[3]);

    printf("avih_chunk_header.four_cc: %c%c%c%c\n",
           m_avih_chunk_header.four_cc.chr[0],
           m_avih_chunk_header.four_cc.chr[1],
           m_avih_chunk_header.four_cc.chr[2],
           m_avih_chunk_header.four_cc.chr[3]);

    printf("avih_chunk_header.size: 0x%x\n\n", m_avih_chunk_header.size);

    printf("main_avih_header.micro_sec_per_frame: %ud\n", m_main_avih_header.micro_sec_per_frame);
    printf("main_avih_header.max_bytes_per_sec: %ud\n", m_main_avih_header.max_bytes_per_sec);
    printf("main_avih_header.padding_granularity: %ud\n", m_main_avih_header.padding_granularity);
    printf("main_avih_header.flags: 0x%x\n", m_main_avih_header.flags);
    printf("main_avih_header.total_frames: %ud\n", m_main_avih_header.total_frames);
    printf("main_avih_header.initial_frames: %ud\n", m_main_avih_header.initial_frames);
    printf("main_avih_header.streams: %ud\n", m_main_avih_header.streams);
    printf("main_avih_header.suggested_buffer_size: 0x%x\n", m_main_avih_header.suggested_buffer_size);
    printf("main_avih_header.width: %d\n", m_main_avih_header.width);
    printf("main_avih_header.height: %d\n\n", m_main_avih_header.height);

    printf("strl_list_header.list: %c%c%c%c\n",
           m_strl_list_header.list.chr[0],
           m_strl_list_header.list.chr[1],
           m_strl_list_header.list.chr[2],
           m_strl_list_header.list.chr[3]);

    printf("strl_list_header.size: 0x%x\n", m_strl_list_header.size);

    printf("strl_list_header.four_cc: %c%c%c%c\n\n",
           m_strl_list_header.four_cc.chr[0],
           m_strl_list_header.four_cc.chr[1],
           m_strl_list_header.four_cc.chr[2],
           m_strl_list_header.four_cc.chr[3]);

    printf("strh_chunk_header.four_cc: %c%c%c%c\n",
           m_strh_chunk_header.four_cc.chr[0],
           m_strh_chunk_header.four_cc.chr[1],
           m_strh_chunk_header.four_cc.chr[2],
           m_strh_chunk_header.four_cc.chr[3]);

    printf("strh_chunk_header.size: 0x%x\n", m_strh_chunk_header.size);

    printf("vids_stream_header.four_cc: %c%c%c%c\n",
           m_vids_stream_header.handler.chr[0],
           m_vids_stream_header.handler.chr[1],
           m_vids_stream_header.handler.chr[2],
           m_vids_stream_header.handler.chr[3]);

    printf("vids_stream_header.flags: 0x%x\n", m_vids_stream_header.flags);
    printf("vids_stream_header.priority: 0x%x\n", m_vids_stream_header.priority);
    printf("vids_stream_header.language: 0x%x\n", m_vids_stream_header.language);
    printf("vids_stream_header.initial_frames: 0x%x\n", m_vids_stream_header.initial_frames);
    printf("vids_stream_header.scale: 0x%x\n", m_vids_stream_header.scale);
    printf("vids_stream_header.rate: %ud\n", m_vids_stream_header.rate);
    printf("vids_stream_header.start: 0x%x\n", m_vids_stream_header.start);
    printf("vids_stream_header.length: 0x%x\n", m_vids_stream_header.length);
    printf("vids_stream_header.suggested_buffer_size: 0x%x\n", m_vids_stream_header.suggested_buffer_size);
    printf("vids_stream_header.quality: 0x%x\n", m_vids_stream_header.quality);
    printf("vids_stream_header.sample_size: 0x%x\n", m_vids_stream_header.sample_size);
    printf("vids_stream_header.frame.left: 0x%x\n", m_vids_stream_header.frame.left);
    printf("vids_stream_header.frame.top: 0x%x\n", m_vids_stream_header.frame.top);
    printf("vids_stream_header.frame.right: 0x%x\n", m_vids_stream_header.frame.right);
    printf("vids_stream_header.frame.bottom: 0x%x\n\n", m_vids_stream_header.frame.bottom);

    printf("strf_chunk_header.four_cc: %c%c%c%c\n",
           m_strf_chunk_header.four_cc.chr[0],
           m_strf_chunk_header.four_cc.chr[1],
           m_strf_chunk_header.four_cc.chr[2],
           m_strf_chunk_header.four_cc.chr[3]);
    printf("strf_chunk_header.size: 0x%x\n\n", m_strf_chunk_header.size);

    printf("bitmap_info_header.size: %ud\n", m_bitmap_info_header.size);
    printf("bitmap_info_header.width: %d\n", m_bitmap_info_header.width);
    printf("bitmap_info_header.height: %d\n", m_bitmap_info_header.height);
    printf("bitmap_info_header.planes: %d\n", m_bitmap_info_header.planes);
    printf("bitmap_info_header.bit_count: %d\n", m_bitmap_info_header.bit_count);
    printf("bitmap_info_header.compression.u32: 0x%x\n", m_bitmap_info_header.compression.u32);
    printf("bitmap_info_header.size_image: %ud\n", m_bitmap_info_header.size_image);
    printf("bitmap_info_header.x_pels_per_meter: %ud\n", m_bitmap_info_header.x_pels_per_meter);
    printf("bitmap_info_header.y_pels_per_meter: %ud\n", m_bitmap_info_header.y_pels_per_meter);
    printf("bitmap_info_header.clr_used: %ud\n", m_bitmap_info_header.clr_used);
    printf("bitmap_info_header.clr_important: %ud\n\n", m_bitmap_info_header.clr_important);

    printf("movi_list_header.list: %c%c%c%c\n",
           m_movi_list_header.list.chr[0],
           m_movi_list_header.list.chr[1],
           m_movi_list_header.list.chr[2],
           m_movi_list_header.list.chr[3]);

    printf("movi_list_header.size: 0x%x\n", m_movi_list_header.size);

    printf("movi_list_header.four_cc: %c%c%c%c\n\n",
           m_movi_list_header.four_cc.chr[0],
           m_movi_list_header.four_cc.chr[1],
           m_movi_list_header.four_cc.chr[2],
           m_movi_list_header.four_cc.chr[3]);

    printf("_00db_chunk_header.four_cc: %c%c%c%c\n",
           m_00db_chunk_header.four_cc.chr[0],
           m_00db_chunk_header.four_cc.chr[1],
           m_00db_chunk_header.four_cc.chr[2],
           m_00db_chunk_header.four_cc.chr[3]);
    printf("_00db_chunk_header.size: 0x%x\n\n", m_00db_chunk_header.size);

    printf("idx1_chunk_header.four_cc: %c%c%c%c\n",
           m_idx1_chunk_header.four_cc.chr[0],
           m_idx1_chunk_header.four_cc.chr[1],
           m_idx1_chunk_header.four_cc.chr[2],
           m_idx1_chunk_header.four_cc.chr[3]);
    printf("idx1_chunk_header.size: 0x%x\n\n", m_idx1_chunk_header.size);

    printf("avi_index_entry.chunk_id: %c%c%c%c\n",
           m_avi_index_entry.chunk_id.chr[0],
           m_avi_index_entry.chunk_id.chr[1],
           m_avi_index_entry.chunk_id.chr[2],
           m_avi_index_entry.chunk_id.chr[3]);
    printf("avi_index_entry.flags: 0x%x\n", m_avi_index_entry.flags);
    printf("avi_index_entry.offset: 0x%x\n", m_avi_index_entry.offset);
    printf("avi_index_entry.size: %ud\n\n", m_avi_index_entry.size);

    return 0;
}

