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
//
// avi_writer.cpp : Defines the exported functions for the DLL application.
//

#include "pipp_video_write.h"  // AVI write class
#include "pipp_avi_write.h"  // AVI write class
#include "pipp_avi_write_dib.h"  // AVI write class
#include "avi_writer.h"
#include <cstdlib>
#include <new>
//#include <mutex>
#include <glib.h>

using namespace std;

#define MAX_CONCURRENT_AVI_FILES 16


// Table of pointer to instances of the AVI write (DIB codec) class
c_pipp_video_write *avi_output_files[MAX_CONCURRENT_AVI_FILES] = {0};

// Mutex to protect the above table
//std::mutex lut_mutex;
GMutex lut_mutex;

// ------------------------------------------
// Create a new AVI file
// ------------------------------------------
int32_t avi_file_create(
    char *filename,
    int32_t  width,
    int32_t  height,
    int32_t  input_format,
	int32_t  codec,
    int32_t  fps)
{
    int ret = -1;

	// Check we support the input format
	int32_t colour;
	if (input_format == AVI_WRITER_INPUT_FORMAT_MONOCHROME) {
		colour = 0;
	} else if (input_format == AVI_WRITER_INPUT_FORMAT_COLOUR) {
		colour = 1;
	} else {
		// Unsupported input format
		return -1;
	}

	// Check we support the requested CODEC
	//if (codec != AVI_WRITER_CODEC_DIB && codec != AVI_WRITER_CODEC_UT_VIDEO) {
    if (codec != AVI_WRITER_CODEC_DIB) {
		// Unsupported codec
		return -1;
	}

	// Check Frames Per Second value
	if (fps <= 0) {
		return -1;
	}

    // Allocate a file ID from the table - start
	int count;
	// Lock mutex to protect table
    //lut_mutex.lock();
	g_mutex_lock(&lut_mutex);
    for (count = 0; count < MAX_CONCURRENT_AVI_FILES; count++) {
        if (avi_output_files[count] == NULL) {
            // This ID slot is free
			// Quickly reserve this slot while under mutex protection;
			avi_output_files[count] = (c_pipp_avi_write_dib *)0xDEADBEEF;

			// Break out of loop
            break;
        }
    }

	// Unlock the mutex so that other threads can access the table
//	lut_mutex.unlock();
	g_mutex_unlock(&lut_mutex);

	// Allocate a file ID from the table - done

	if (count == MAX_CONCURRENT_AVI_FILES) {
		// Could not allocate slot
		return -1;
	}

	// Create a new instance of the required AVI video class
    if (codec == AVI_WRITER_CODEC_DIB) {
        avi_output_files[count] = new(std::nothrow) c_pipp_avi_write_dib;
    }// else if (codec == AVI_WRITER_CODEC_UT_VIDEO) {
//        avi_output_files[count] = new(std::nothrow) c_pipp_avi_write_utvideo;
//    }

	// Check if the new instance was actually created
	if (avi_output_files[count] == NULL) {
		// The new operator failed
		// No cleanup required
		return -1;
	}

	// Actually create the AVI file ready for use
	ret = avi_output_files[count]->create(
		filename,
		width,
		height,
		colour,
		fps,
        1,
		0,
		0);

	if (ret == 0) {
		// The file was correctly created, return FILE ID
		return count;
	} else {
		// The file could not be created
		delete avi_output_files[count];  // Delete class instance
		avi_output_files[count] = NULL;  // Free up slot in table
		return ret;
	}
}


// ------------------------------------------
// Write frame to AVI file
// ------------------------------------------
int32_t avi_file_write_frame(
            int32_t file_id,
            uint8_t *data)
{
	// Check file ID is valid
    if (avi_output_files[file_id] == NULL || file_id >= MAX_CONCURRENT_AVI_FILES) {
        return -1;
    }

    return avi_output_files[file_id]->write_frame(
        data,
        0,
        1);
}


// ------------------------------------------
// Close the AVI file
// ------------------------------------------
int32_t avi_file_close(int32_t file_id)
{
	// Check file ID is valid
    if (avi_output_files[file_id] == NULL || file_id >= MAX_CONCURRENT_AVI_FILES) {
        return -1;
    }

	int32_t ret = avi_output_files[file_id]->close();
	delete avi_output_files[file_id];  // Free c_pipp_avi_write_dib class instance
	avi_output_files[file_id] = NULL;

    return ret;
}
