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


#ifndef PIPP_BUFFER_H
#define PIPP_BUFFER_H

#include <cstdint>
#include <cstdint>
#include <cstddef>


//#define DEBUG


class c_pipp_buffer {
    // ------------------------------------------
    // Private definitions
    // ------------------------------------------
    private:
        uint32_t m_buffer_size;
        uint8_t *mp_buffer;

    // ------------------------------------------
    // Public definitions
    // ------------------------------------------

    public:
        // ------------------------------------------
        // Constructor
        // ------------------------------------------
        c_pipp_buffer() : 
            m_buffer_size(0),
            mp_buffer(nullptr) {};


        // ------------------------------------------
        // Destructor
        // ------------------------------------------
        ~c_pipp_buffer() {
            delete [] mp_buffer;
        };


        // ------------------------------------------
        // Get buffer ptr
        // ------------------------------------------
        inline uint8_t *get_buffer_ptr() {
            return mp_buffer;
        };


        // ------------------------------------------
        // Member function to set a new  buffer
        // ------------------------------------------
        // Caller has the responsibility to esure the old buffer is freed
        // and the new buffer is large enough
        uint8_t *set_buffer_ptr(uint8_t *new_ptr) {
            uint8_t *ret = mp_buffer;
            mp_buffer = new_ptr;
            return ret;
        }


        // ------------------------------------------
        // Member function to get a new  buffer
        // ------------------------------------------
        uint8_t *get_buffer(
            uint32_t size);


        // ------------------------------------------
        // Member function to get a new zeroed buffer
        // ------------------------------------------
        uint8_t *get_zero_buffer(
            uint32_t size);


        // ------------------------------------------
        // Member function to delete buffer
        // ------------------------------------------
        void delete_buffer();

};


#endif  // PIPP_BUFFER_H
