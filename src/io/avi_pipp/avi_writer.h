// ---------------------------------------------------------------------
// Copyright (C) 2015 Chris Garry
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>
// ---------------------------------------------------------------------

#ifndef SRC_IO_AVI_WRITER_H_
#define SRC_IO_AVI_WRITER_H_


#include <stdint.h>

#define AVI_WRITER_INPUT_FORMAT_MONOCHROME  0
#define AVI_WRITER_INPUT_FORMAT_COLOUR      1
#define AVI_WRITER_CODEC_DIB       0
#define AVI_WRITER_CODEC_UT_VIDEO  1

#ifdef __cplusplus
extern "C" {
#endif

    // Create a new AVI file
int32_t avi_file_create(
        char *filename,
        int32_t  width,          // Frame width in pixels
        int32_t  height,         // Frame height in pixels
        int32_t  input_format,   // 0: Monochrome LLL, 1: Colour BGRBGRBGR
        int32_t  codec,          // 0: DIB, 1: Ut Video
        int32_t  fps);           // Frames per second value


    // Write frame to video file
int32_t avi_file_write_frame(int32_t file_id, uint8_t *data);

    // Close the AVI file
int32_t avi_file_close(int32_t file_id);


#ifdef __cplusplus
}
#endif

#endif /* SRC_IO_AVI_WRITER_H_ */
