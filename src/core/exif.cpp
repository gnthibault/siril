/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2019 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __cplusplus
extern "C" {
#endif
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <glib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __cplusplus
}
#endif

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#if defined(_WIN32) && defined(EXV_UNICODE_PATH)
  #define WIDEN(s) pugi::as_wide(s)
#else
#define WIDEN(s) (s)
#endif

#include <exiv2/exiv2.hpp>

#include "exif.h"

/** code from darktable */
/**
 * Get the largest possible thumbnail from the image
 */
int siril_get_thumbnail_exiv(const char *path, uint8_t **buffer, size_t *size,
		char **mime_type) {
	try {
		std::unique_ptr<Exiv2::Image> image(
				Exiv2::ImageFactory::open(WIDEN(path)));
		assert(image.get() != 0);
		image->readMetadata();

		// Get a list of preview images available in the image. The list is sorted
		// by the preview image pixel size, starting with the smallest preview.
		Exiv2::PreviewManager loader(*image);
		Exiv2::PreviewPropertiesList list = loader.getPreviewProperties();
		if (list.empty()) {
			std::cerr << "[exiv2] couldn't find thumbnail for " << path << std::endl;
			return 1;
		}

		// Select the largest one
		// FIXME: We could probably select a smaller thumbnail to match the mip size
		//        we actually want to create. Is it really much faster though?
		Exiv2::PreviewProperties selected = list.back();

		// Get the selected preview image
		Exiv2::PreviewImage preview = loader.getPreviewImage(selected);
		const unsigned char *tmp = preview.pData();
		size_t _size = preview.size();

		*size = _size;
		*mime_type = strdup(preview.mimeType().c_str());
		*buffer = (uint8_t*) malloc(_size);
		if (!*buffer) {
			std::cerr << "[exiv2] couldn't allocate memory for thumbnail for "
					<< path << std::endl;
			return 1;
		}
		//std::cerr << "[exiv2] "<< path << ": found thumbnail "<< preview.width() << "x" << preview.height() << std::endl;
		memcpy(*buffer, tmp, _size);

		return 0;
	} catch (Exiv2::AnyError &e) {
		std::string s(e.what());
		std::cerr << "[exiv2]: " << s << std::endl;
		return 1;
	}
}
