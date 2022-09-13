/*
 * LUT.h
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2011 Jan Rinze Peterzon (janrinze@gmail.com)
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

/*
 *  Declaration of flexible Lookup Tables
 *
 *  Usage:
 *
 *      LUT<type> name (size);
 *      LUT<type> name (size, flags);
 *
 *      creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *          LUT<float> my_lut (10);
 *          float value = my_lut[3];
 *          float value = my_lut[2.5]; // this will interpolate
 *
 *      when using a float type index it will interpolate the lookup values
 *
 *      extra setting in flags: (clipping is set by default)
 *      LUT_CLIP_ABOVE
 *      LUT_CLIP_BELOW
 *
 *      example:
 *          LUT<float> my_lut (10,LUT_CLIP_BELOW);
 *          float value = my_lut[22.5];  // this will extrapolate
 *          float value = my_lut[-22.5]; // this will not extrapolate
 *
 *          LUT<float> my_lut (10,0); // this will extrapolate on either side
 *
 *      shotcuts:
 *
 *          LUTf stands for LUT<float>
 *          LUTi stands for LUT<int>
 *          LUTu stands for LUT<unsigned int>
 *          LUTd stands for LUT<double>
 *          LUTuc stands for LUT<unsigned char>
 */


#include <algorithm>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <vector>

#ifndef NDEBUG
#include <fstream>
#endif

#include "opthelper.h"
#include "rt_math.h"

// Bit representations of flags
enum {
    LUT_CLIP_OFF,   // LUT does not clip input values
    LUT_CLIP_BELOW, // LUT clips input values at lower bound
    LUT_CLIP_ABOVE  // LUT clips input values at upper bound
};

template<typename T>
class LUT;

using LUTf = LUT<float>;
using LUTi = LUT<int32_t>;
using LUTu = LUT<uint32_t>;
using LUTd = LUT<double>;
using LUTuc = LUT<uint8_t>;

template<typename T>
class LUT
{
protected:
    // list of variables ordered to improve cache speed
    int maxs;
    float maxsf;
    T * data;
    unsigned int clip;
    unsigned int size;
    unsigned int upperBound;  // always equals size-1, parameter created for performance reason
private:
    unsigned int owner;
public:
    /// convenience flag! If one doesn't want to delete the buffer but want to flag it to be recomputed...
    /// The user have to handle it itself, even if some method can (re)initialize it
    bool dirty;

    explicit LUT(int s, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE, bool initZero = false)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);
#endif
        dirty = true;
        clip = flags;
        // Add a few extra elements so [](vfloat) won't access out-of-bounds memory.
        // The routine would still produce the right answer, but might cause issues
        // with address/heap checking programs.
        data = new T[s + 3];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
        if (initZero) {
            clear();
        }
    }

    explicit LUT(const std::vector<T>& input, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE) :
        maxs(input.size() - 2),
        maxsf(maxs),
        data(new T[input.size() + 3]), // Add a few extra elements so [](vfloat) won't access out-of-bounds memory.
        clip(flags),
        size(input.size()),
        upperBound(size - 1),
        owner(1),
        dirty(true)
    {
#ifndef NDEBUG

        if (input.empty()) {
            printf("s=0!\n");
        }

        assert(!input.empty());
#endif
        std::copy_n(input.begin(), input.size(), data);
    }

    void operator ()(int s, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE, bool initZero = false)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);
#endif

        if (owner && data) {
            delete[] data;
        }

        dirty = true; // Assumption!
        clip = flags;
        // See comment in constructor.
        data = new T[s + 3];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
        if (initZero) {
            clear();
        }

    }

    LUT()
    {
        data = nullptr;
        reset();
    }

    ~LUT()
    {
        if (owner) {
            delete[] data;
#ifndef NDEBUG
            data = (T*)0xBAADF00D;
#endif
        }
    }

    explicit LUT(const LUT&) = delete;

    void setClip(int flags)
    {
        clip = flags;
    }

    int getClip() const {
        return clip;
    }

    /** @brief Get the number of element in the LUT (i.e. dimension of the array)
     *  For a LUT(500), it will return 500
     *  @return number of element in the array
     */
    unsigned int getSize() const
    {
        return size;
    }

    /** @brief Get the highest value possible (i.e. dimension of the array)
     *  For a LUT(500), it will return 499, because 500 elements, starting from 0, goes up to 499
     *  @return number of element in the array
     */
    unsigned int getUpperBound() const
    {
        return size > 0 ? upperBound : 0;
    }

    LUT<T>& operator=(const LUT<T>& rhs)
    {
        if (this != &rhs) {
            if (rhs.size > this->size) {
                delete [] this->data;
                this->data = nullptr;
            }

            if (this->data == nullptr) {
                // See comment in constructor.
                this->data = new T[rhs.size + 3];
            }

            this->clip = rhs.clip;
            this->owner = 1;
            memcpy(this->data, rhs.data, rhs.size * sizeof(T));
            this->size = rhs.size;
            this->upperBound = rhs.upperBound;
            this->maxs = this->size - 2;
            this->maxsf = (float)this->maxs;
        }

        return *this;
    }

    // handy to sum up per thread histograms. //#pragma omp simd speeds up the loop by about factor 3 for LUTu (uint32_t).
    LUT<T>& operator+=(const LUT<T>& rhs)
    {
        if (rhs.size == this->size) {
#ifdef _OPENMP
            //#pragma omp simd
#endif

            for(unsigned int i = 0; i < this->size; i++) {
                data[i] += rhs.data[i];
            }
        }

        return *this;
    }

    // multiply all elements of LUT<float> with a constant float value
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    LUT<float>& operator*=(float factor)
    {
#ifdef _OPENMP
        //#pragma omp simd
#endif

        for(unsigned int i = 0; i < this->size; i++) {
            data[i] *= factor;
        }

        return *this;
    }

    // divide all elements of LUT<float> by a constant float value
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    LUT<float>& operator/=(float divisor)
    {
#ifdef _OPENMP
        //#pragma omp simd
#endif

        for(unsigned int i = 0; i < this->size; i++) {
            data[i] /= divisor;
        }

        return *this;
    }


    // use with integer indices
    T& operator[](int index) const
    {
        return data[ LIM<int>(index, 0, upperBound) ];
    }


    // use with float indices
    template<typename U = T, typename V, typename = typename std::enable_if<std::is_floating_point<V>::value && std::is_same<U, float>::value>::type>
    T operator[](V index) const
    {
        int idx = (int)index;  // don't use floor! The difference in negative space is no problems here

        if (index < 0.f) {
            if (clip & LUT_CLIP_BELOW) {
                return data[0];
            }

            idx = 0;
        } else if (index > maxsf) {
            if (clip & LUT_CLIP_ABOVE) {
                return data[upperBound];
            }

            idx = maxs;
        }

        float diff = index - (float) idx;
        T p1 = data[idx];
        T p2 = data[idx + 1] - p1;
        return (p1 + p2 * diff);
    }

    // Return the value for "index" that is in the [0-1] range.
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    T getVal01(float index) const
    {
        index *= (float)upperBound;
        int idx = (int)index;  // don't use floor! The difference in negative space is no problems here

        if (index < 0.f) {
            if (clip & LUT_CLIP_BELOW) {
                return data[0];
            }

            idx = 0;
        } else if (index > maxsf) {
            if (clip & LUT_CLIP_ABOVE) {
                return data[upperBound];
            }

            idx = maxs;
        }

        float diff = index - (float) idx;
        T p1 = data[idx];
        T p2 = data[idx + 1] - p1;
        return (p1 + p2 * diff);
    }

    operator bool() const // FIXME: Should be explicit
    {
        return size > 0;
    }

    void clear()
    {
        if (data && size) {
            memset(data, 0, size * sizeof(T));
        }
    }

    void reset()
    {
        if (data) {
            delete[] data;
        }

        dirty = true;
        data = nullptr;
        owner = 1;
        size = 0;
        upperBound = 0;
        maxs = 0;
        maxsf = 0.f;
        clip = 0;
    }

    // create an identity LUT (LUT(x) = x) or a scaled identity LUT (LUT(x) = x / divisor)
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    void makeIdentity(float divisor = 1.f)
    {
        if(divisor == 1.f) {
            for(unsigned int i = 0; i < size; i++) {
                data[i] = i;
            }
        } else {
            for(unsigned int i = 0; i < size; i++) {
                data[i] = i / divisor;
            }
        }
    }

    // compress a LUT<uint32_t> with size y into a LUT<uint32_t> with size x (y>x)
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void compressTo(LUT<T> &dest, unsigned int numVals = 0) const
    {
        numVals = numVals == 0 ? size : numVals;
        numVals = std::min(numVals, size);
        float divisor = numVals - 1;
        float mult = (dest.size - 1) / divisor;

        for (unsigned int i = 0; i < numVals; i++) {
            int hi = (int)(mult * i);
            dest.data[hi] += this->data[i] ;
        }
    }

    // compress a LUT<uint32_t> with size y into a LUT<uint32_t> with size x (y>x) by using the passThrough LUT to calculate indexes
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void compressTo(LUT<T> &dest, unsigned int numVals, const LUT<float> &passThrough) const
    {
        if(passThrough) {
            numVals = std::min(numVals, size);
            numVals = std::min(numVals, passThrough.getSize());
            float mult = dest.size - 1;

            for (unsigned int i = 0; i < numVals; i++) {
                int hi = (int)(mult * passThrough[i]);
                dest[hi] += this->data[i] ;
            }
        }
    }

    // compute sum and average of a LUT<uint32_t>
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void getSumAndAverage(float &sum, float &avg) const
    {
        sum = 0.f;
        avg = 0.f;
        int i = 0;


        for (; i < static_cast<int>(size); i++) {
            T val = data[i];
            sum += val;
            avg += i * val;
        }

        avg /= sum;
    }


    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    void makeConstant(float value, unsigned int numVals = 0)
    {
        numVals = numVals == 0 ? size : numVals;
        numVals = std::min(numVals, size);

        for(unsigned int i = 0; i < numVals; i++) {
            data[i] = value;
        }
    }

    // share the buffer with another LUT, handy for same data but different clip flags
    void share(const LUT<T> &source, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE)
    {
        if (owner && data) {
            delete[] data;
        }

        dirty = false;  // Assumption
        clip = flags;
        data = source.data;
        owner = 0;
        size = source.getSize();
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
    }


};
