// Implementation of the technique described in the blog post "Optimizing binary search"
// by David Geier (visit http://geidav.wordpress.com)
#define NDEBUG 
#include <functional>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <random>
#include <vector>
#include <cstdint>

// mapping function: specialized for 64-bit signed/unsigned
// integers and 64-bit floating point values.
// the mapping functions are used to create the LUT, because signed
// integers and floats are not comparable using bit-wise comparison.



template<class T> uint64_t MapValue(T val)
{
    // default implementation results in a compile-time error
    static_assert(sizeof(T) == 0, "type unavailable: only 64-bit signed/unsigned int and float supported");
    return 0;
}

template<> uint64_t MapValue<uint64_t>(uint64_t val) // 64-bit unsigned int
{
    return val; // no bit-twiddling required, just forward
}

template<> uint64_t MapValue<int64_t>(int64_t val) // 64-bit signed int
{
    return (uint64_t)val^0x8000000000000000; // flip sign bit
}

// taken from Michael Herf's article on Stereopsis


// LUT optimized binary search implementation for 64-bit POD types
template<class T, size_t LUT_BITS> class SearchPod64
{
public:
    SearchPod64(const std::vector<T> &vals) :
        Vals(vals)
    {
        static_assert(LUT_BITS > 0 && LUT_BITS < 64, "invalid binary search LUT size");
        InitLut();
    }

    ssize_t StdBinarySearch(T key) const
    {
        const auto iter = std::lower_bound(Vals.begin(), Vals.end(), key);
        return (iter != Vals.end() && *iter == key ? std::distance(Vals.begin(), iter) : -1);
    }

    ssize_t MyBinarySearch(T key) const
    {
        return BinarySearch(0, (ssize_t)Vals.size()-1, key);
    }

    ssize_t LutBinarySearch(T key) const
    {
        if (key <= Vals.front()) return 0;
        if (key > Vals.back()) return Vals.size();

        const auto mappedKey = MapValue<T>(key);
        const auto lutIdx = mappedKey>>(64-LUT_BITS);
        const auto start = Lut[lutIdx];

        // interval end of i-th LUT entry = interval start of (i+1)th LUT entry - 1.
        // however, all LUT remaining LUT entries map to the last valid interval
        // start => interval end < interval start => just use number of values
        const auto end = (lutIdx+1 >= LutEnd ? Vals.size()-1 : Lut[lutIdx+1]-1);
        return BinarySearch(start, end, key);
    }
    ssize_t UnsafeLutBinarySearch(T key) const
    {
        const auto mappedKey = MapValue<T>(key);
        const auto lutIdx = mappedKey>>(64-LUT_BITS);
        const auto start = Lut[lutIdx];

        // interval end of i-th LUT entry = interval start of (i+1)th LUT entry - 1.
        // however, all LUT remaining LUT entries map to the last valid interval
        // start => interval end < interval start => just use number of values
        const auto end = (lutIdx+1 >= LutEnd ? Vals.size()-1 : Lut[lutIdx+1]-1);
        return BinarySearch(start, end, key);
    }
private:
    ssize_t BinarySearch(ssize_t left, ssize_t right, T key) const
    {
        if (key > Vals[right]) return right + 1; 
        
        /*
        size_t __len = right-left;
        size_t __first = left;
        while (__len > 0)
        {
            size_t __half = __len >> 1;
            size_t __middle = __first+__half;
            if (Vals[__middle] < key)
            {
                __first = __middle;
                ++__first;
                __len = __len - __half - 1;
            }
            else
                __len = __half;
        }
        return __first;
        */
        while (left < right)
        {
            const auto mid = left+((right-left)>>1); // avoids overflow
            assert(mid < right); // interval must be reduced in each iteration
            const auto valMid = Vals[mid];

            // no early exit so that always the occurence
            // of the key with the lowest index is found
            if (valMid < key)
                left = mid+1;
            else
                right = mid;
        }

        assert(left == right);
        return left; 
    }

    void InitLut()
    {
        // fill look-up-table
        Lut.resize((1<<LUT_BITS)+1); // one additional element to avoid condition in interval end computation

        size_t thresh = 0;
        size_t last = 0;

        for (ssize_t i=0; i<(ssize_t)Vals.size()-1; i++)
        {
            const uint64_t mappedNextVal = MapValue<T>(Vals[i+1]);
            const uint64_t nextThresh = mappedNextVal>>(64-LUT_BITS);
            Lut[thresh] = last;

            if (nextThresh > thresh)
            {
                last = i+1;
                for (size_t j=thresh+1; j<=nextThresh; j++)
                    Lut[j] = last;
            }

            thresh = nextThresh;
        }

        // set remaining thresholds that couldn't be found
        for (size_t i=thresh; i<Lut.size()-1; i++)
            Lut[i] = last;

        // remember last end threshold index, because the interval
        // end of all values mapping to an entry bigger than
        // that have to be handled differently
        LutEnd = thresh;

        /*
        for (auto v : Lut)
            std::cout << v.first << ", " << v.second << std::endl;
        */
    }

private:
    std::vector<size_t>    Lut;
    const std::vector<T> & Vals;
    size_t                 LutEnd;
};




void int64BinarySearch(uint64_t* data, long int N, uint64_t* keys, long int M, int32_t* result, uint64_t mult)
{
    {
        std::vector<uint64_t> keyVec(M);
        for (int i=0; i<M; i++)
        {
        keyVec[i] = keys[i]*mult;
        }
        SearchPod64<uint64_t, 24> s(keyVec);
        
        #pragma omp parallel for shared(s,data,result,mult)  num_threads(8) 
        for (int i =0; i<N; i++)
        {
            auto ids = s.LutBinarySearch(data[i]*mult);    
            result[i] = (int32_t) ids;       
        }
    }

    return; 
}

/*
void uint64BinarySearch(uint64_t* data, long int N, uint64_t* keys, long int M, int32_t* result, uint64_t multiplyer){
    std::vector<int64_t> keyVec(M);
    for (int i=0; i<M; i++)
    {
    keyVec[i] = keys[i] * multiplyer;
    }
    SearchPod64<int64_t, 24> s(keyVec);
    
    //
    for (int i =0; i<N; i++)
    {
        auto ids = s.LutBinarySearch(data[i] * multiplyer);    
        result[i] = (int32_t) ids;       
    }
    return; 
}


*/
