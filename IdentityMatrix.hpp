/**
 * @file mtl_test.cpp
 * Implimentation file for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include<iostream>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL

class IdentityMatrix{
    public:
    IdentityMatrix(unsigned int size):n(size){}
    

    unsigned int get_size() const
    {
        return n;
    }

    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn &v, VectorOut &w, Assign) const
    {
        Assign::apply(w, v);
    }

    template <typename Vector>
    mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>    
    operator*(const Vector &v) const
    {
        return {*this, v};
    }

    private: 
    unsigned int n;
};

inline std::size_t size(const IdentityMatrix A)
{
    return A.get_size()*A.get_size();
}
inline std::size_t num_rows(const IdentityMatrix A)
{
    return A.get_size();
}
inline std::size_t num_cols(const IdentityMatrix A)
{
    return A.get_size();
}


namespace mtl{
    namespace ashape{
        template<>
        struct ashape_aux<IdentityMatrix>{
            typedef nonscal type;
        };
    }//

    template<>
    struct Collection<IdentityMatrix>{
        typedef double value_type;
        typedef unsigned size_type;
    };
};