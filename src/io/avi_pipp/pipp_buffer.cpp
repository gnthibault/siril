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


#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include "pipp_buffer.h"


using namespace std;


// ------------------------------------------
// Member function to get a new buffer
// ------------------------------------------
uint8_t *c_pipp_buffer::get_buffer(
    uint32_t size)
{
    if (size > m_buffer_size) {
        // The current buffer is not big enough
        // Delete the current one and generate a new one that is large enough
        delete [] mp_buffer;

        // Calculate buffer size
        m_buffer_size = 128;
        while (m_buffer_size < size) {
            m_buffer_size *= 2;
        }

        // Create the buffer
        try {
            mp_buffer = new uint8_t[m_buffer_size];
        } catch(...) {
            cout << "FATAL ERROR: memory allocation (";
            cout << dec << (size * sizeof(m_buffer_size));
            cout << ") failed in c_pipp_buffer::get_buffer()" << endl;
            exit(-1);
        }
    }

    // Return pointer to buffer
    return mp_buffer;
}


// ------------------------------------------
// Member function to get a new zeroed buffer
// ------------------------------------------
uint8_t *c_pipp_buffer::get_zero_buffer(
    uint32_t size)
{
    // Call get_buffer() method to create buffer as usual
    get_buffer(size);

    // Zero the buffer
    memset(mp_buffer, 0, m_buffer_size);

    // Return pointer to buffer
    return mp_buffer;
}


// ------------------------------------------
// Member function to delete buffer
// ------------------------------------------
void c_pipp_buffer::delete_buffer() {
    delete [] mp_buffer;
    mp_buffer = nullptr;
    m_buffer_size = 0;
}
