
#include "div_by_s_i.hpp"
namespace cantor {

    template<typename FieldT>   
    std::vector<FieldT> additive_FFT_hc(const std::vector<FieldT> &poly_coeffs, const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        // size_t diff_dim = m-log_poly_terms;
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull<<m));

        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, 192, or 256 for using the cantor basis");

        unsigned index = 1<<log_poly_terms;
        while(index < n){
            std::copy(poly_coeffs.begin(), poly_coeffs.end(), g.begin() + index);
            index += (1<<log_poly_terms);
        }

        const unsigned* nz_S;
        size_t S_index=log_poly_terms, t; 
        size_t input_size = 1<<log_poly_terms;
        size_t n_modules = 1<< (m - log_poly_terms);
        size_t half_input_size = input_size >> 1;
        size_t offset = 0;
        size_t shift_bit;

        
        size_t m_prime = log_poly_terms;

        if(m_prime>9){
            for (int r = 0; r < m_prime-9; ++r)
            {   
                --S_index; 
                offset = 0;
                nz_S = libiop::s_i[S_index]; t=libiop::n_terms[S_index];
                shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);

                for (size_t module = 0; module < n_modules; ++module){
                    size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                    FieldT mult_factor = FieldT::zero();
                    while(module_shifted){
                        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                        i_256 += 256;
                        module_shifted >>= 8;
                    }
                    size_t offset2 = offset + half_input_size;
                    for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                        FieldT gk = g[k];
                        for (unsigned i = 0; i < t; ++i){
                            g[k - nz_S[i]] += gk;
                        }
                        g[k-half_input_size] += gk * mult_factor;
                    }
                    for (size_t j = 0; j < half_input_size; ++j)
                        g[offset2+j] += g[offset+j];

                    offset += input_size;
                }
                input_size = half_input_size;
                n_modules <<= 1;
                half_input_size = input_size >> 1;
            }
        }
        
        if(m_prime>8){    // r = m-9 (S_8) (input_size = 512)
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }
                DIV_S8(mult_factor,offset);
                offset += input_size;
            }
            input_size = half_input_size;
            n_modules <<= 1;
            S_index = 8; 
        }

        if(m_prime>5){
            offset = 0;
            half_input_size = input_size >> 1;
            for (int r = m-8; r < m-5; ++r)
            {   
                shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
                offset = 0;
                --S_index; 
                nz_S = libiop::s_i[S_index]; t=libiop::n_terms[S_index];
                for (size_t module = 0; module < n_modules; ++module){
                    size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                    FieldT mult_factor = FieldT::zero();
                    while(module_shifted){
                        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                        i_256 += 256;
                        module_shifted >>= 8;
                    }
                    size_t offset2 = offset + half_input_size;
                    for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                        FieldT gk = g[k];
                        for (unsigned i = 0; i < t; ++i){
                            g[k - nz_S[i]] += gk;
                        }
                        g[k-half_input_size] += gk * mult_factor;
                    }
                    for (size_t j = 0; j < half_input_size; ++j)
                        g[offset2+j] += g[offset+j];

                    offset += input_size;
                }
                input_size = half_input_size;
                half_input_size = input_size >> 1;
                n_modules <<= 1;
            }
        }

        if(m_prime>4){     // r = m-5 (S_4) (input_size = 32)
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }
                DIV_S4(mult_factor,offset);
                offset += 32;
            }
            input_size = half_input_size;
            n_modules <<= 1;
        }

        if(m_prime>3){ // r = m-4 (S_3) (input_size = 16)
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                DIV_S3(mult_factor,offset);
                offset += 16;
            }
            input_size = half_input_size;
            n_modules <<= 1;
        }

        if(m_prime>2){ // r = m-3 (S_2) (input_size = 8)
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                DIV_S2(mult_factor,offset);
                offset += 8;
            }
            input_size = half_input_size;
            n_modules <<= 1;
        }

        if(m_prime>1){ // r = m-2 (S_1) (input_size = 4)
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                DIV_S1(mult_factor,offset);
                offset += 4;
            }
            input_size = half_input_size;
            n_modules <<= 1;
        }

        if(m_prime>0){ // r = m-1 (S_0) (input_size = 4)    
            offset = 0;
            half_input_size = input_size >> 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                DIV_S0(mult_factor,offset);
                offset += 2;
            }
        }

        return g;
    }

    template<typename FieldT>   
    std::vector<FieldT> additive_IFFT_hc(const std::vector<FieldT> &evals, const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(evals);
        const size_t n = g.size();
        assert(n == (1ull<<m));

        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) libiop::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, or 256 for using the cantor basis");

        size_t input_size = 1;
        size_t n_modules = n;

        const unsigned* nz_S;
        size_t S_index = 0, t;

        size_t offset = 0;
        size_t half_input_size;
        size_t shift_bit;

        if(m>0){ // r = m-1 (S_0) (input_size = 2)   
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                INV_DIV_S0(mult_factor,offset);
                offset += 2;
            }
        }

        if(m>1){ // r = m-2 (S_1) (input_size = 4)
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                INV_DIV_S1(mult_factor,offset);
                offset += 4;
            }
        }

        if(m>2){ // r = m-3 (S_2) (input_size = 8)
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                INV_DIV_S2(mult_factor,offset);
                offset += 8;
            }
        }

        if(m>3){ // r = m-4 (S_3) (input_size = 16)
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }                
                INV_DIV_S3(mult_factor,offset);
                offset += 16;
            }
        }

        if(m>4){     // r = m-5 (S_4) (input_size = 32)
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){    
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }
                INV_DIV_S4(mult_factor,offset);
                offset += 32;
            }
        }

        if(m>5){
            S_index = 5;
            // half_input_size = input_size >> 1;
            while (S_index >= 5 && S_index < 8)
            {   
                nz_S = libiop::s_i[S_index]; t=libiop::n_terms[S_index];
                offset = 0;
                half_input_size = input_size;
                n_modules >>= 1;
                input_size <<= 1;
                shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
                for (size_t module = 0; module < n_modules; ++module){
                    size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                    FieldT mult_factor = FieldT::zero();
                    while(module_shifted){
                        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                        i_256 += 256;
                        module_shifted >>= 8;
                    }
                    size_t offset2 = offset + half_input_size;
                    for (size_t j = 0; j < half_input_size; ++j)
                        g[offset2+j] += g[offset+j];
                    for (size_t k =  offset2; k < offset2+half_input_size; ++k){
                        FieldT gk = g[k];
                        g[k-half_input_size] += gk * mult_factor;
                        for (size_t i = 0; i < t; ++i){
                            g[k - nz_S[i]] += gk;
                        }
                    }
                    offset += input_size;
                }
                ++S_index; 
            }
        }


        if(m>8){    // r = m-9 (S_8) (input_size = 512)
            offset = 0;
            half_input_size = input_size;
            input_size <<= 1;
            n_modules >>= 1;
            shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
            for (unsigned module = 0; module < n_modules; ++module){
                size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                FieldT mult_factor = FieldT::zero();
                while(module_shifted){
                    mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                    i_256 += 256;
                    module_shifted >>= 8;
                }
                INV_DIV_S8(mult_factor,offset);
                offset += input_size;
            }
        }

        if(m>9){
            S_index = 9;
            while (S_index >= 9 && S_index < m)
            {   
                nz_S = libiop::s_i[S_index]; t=libiop::n_terms[S_index];
                offset = 0;
                half_input_size = input_size;
                n_modules >>= 1;
                input_size <<= 1;
                shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
                for (size_t module = 0; module < n_modules; ++module){
                    size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
                    FieldT mult_factor = FieldT::zero();
                    while(module_shifted){
                        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                        i_256 += 256;
                        module_shifted >>= 8;
                    }
                    size_t offset2 = offset + half_input_size;
                    for (size_t j = 0; j < half_input_size; ++j)
                        g[offset2+j] += g[offset+j];
                    for (size_t k =  offset2; k < offset2+half_input_size; ++k){
                        FieldT gk = g[k];
                        g[k-half_input_size] += gk * mult_factor;
                        for (size_t i = 0; i < t; ++i){
                            g[k - nz_S[i]] += gk;
                        }
                    }
                    offset += input_size;
                }
                ++S_index; 
            }
        }

        return g;
    }


} // namespace cantor
