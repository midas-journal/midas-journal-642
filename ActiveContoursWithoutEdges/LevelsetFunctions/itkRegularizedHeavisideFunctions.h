#ifndef __itkRegularizedHeavisideFunction_h
#define __itkRegularizedHeavisideFunction_h

#include <vcl_cmath.h>
#include <vnl/vnl_math.h>

namespace itk
{
  /** \class HeavisideFunctionBase
  * \brief Implementation of theoretical definitions of the Heaviside and dirac functions
  */
  template< typename TOutput = double, typename TInput = float >
  class HeavisideFunctionBase
  {
  public:
    typedef TInput InputType;
    typedef TOutput OutputType;
    typedef HeavisideFunctionBase Self;

    HeavisideFunctionBase() : m_Epsilon( 1. ),
      m_InvEpsilon( 1. )
    {}
    HeavisideFunctionBase( const OutputType& iEps ) : m_Epsilon( iEps ),
      m_InvEpsilon( 1./iEps )
    {}
    HeavisideFunctionBase( const Self& iOther ) : m_Epsilon( iOther.m_Epsilon ),
      m_InvEpsilon( 1./ iOther.m_Epsilon )
    {}

    virtual ~HeavisideFunctionBase() {}

    void SetEpsilon( const OutputType& iEps )
    {
      m_Epsilon = iEps;

      if ( iEps > 0 )
      {
        m_InvEpsilon = 1./iEps;
      }
      else
      {
        std::cerr << "ERROR: Epsilon needs to be greater than 0" << std::endl;
      }
    }

    OutputType GetEpsilon() const
    {
      return m_Epsilon;
    }

    /** \brief Heaviside Function */
    virtual OutputType Heaviside( const InputType& iX ) const
    {
      return ( iX >= 0. ) ? 1. : 0.;
    }

    /** \brief Dirac Function */
    virtual OutputType Dirac( const InputType& iX ) const
    {
      return ( iX == 0. ) ? 1. : 0.;
    }

  protected:
    OutputType m_Epsilon;
    OutputType m_InvEpsilon;
  };

  template< typename TOutput = double, typename TInput = float >
  class AtanHeavisideFunction : public HeavisideFunctionBase< TOutput, TInput >
  {
  public:
    typedef TInput InputType;
    typedef TOutput OutputType;
    typedef HeavisideFunctionBase< TOutput, TInput > Superclass;
    typedef AtanHeavisideFunction Self;

    AtanHeavisideFunction() : Superclass()
    {}
    AtanHeavisideFunction( const OutputType& iEps ) : Superclass( iEps )
    {}
    ~AtanHeavisideFunction()
    {}

    inline OutputType Heaviside( const InputType& iX ) const
    {
      return 0.5 + ( vnl_math::one_over_pi * vcl_atan( iX * this->m_InvEpsilon ) );
    }

    inline OutputType Dirac( const InputType& iX ) const
    {
      OutputType t = ( iX * this->m_InvEpsilon );
      return vnl_math::one_over_pi * (1. + t * t );
    }
  };

  template< typename TOutput = double, typename TInput = float >
  class SinHeavisideFunction : public HeavisideFunctionBase< TOutput, TInput >
  {
  public:
    typedef TInput InputType;
    typedef TOutput OutputType;
    typedef HeavisideFunctionBase< TOutput, TInput > Superclass;
    typedef SinHeavisideFunction Self;

    SinHeavisideFunction() : Superclass()
    {}
    SinHeavisideFunction( const OutputType& iEps ) : Superclass( iEps )
    {}
    ~SinHeavisideFunction()
    {}

    inline OutputType Heaviside( const InputType& iX ) const
    {
      if( iX > this->m_Epsilon )
      {
        return 1.;
      }
      else
      {
        if( iX < -this->m_Epsilon )
        {
          return 0.;
        }
        else
        {
          OutputType t = iX * this->m_InvEpsilon;
          return 0.5 + ( t + vcl_sin( vnl_math::pi * iX ) * vnl_math::one_over_pi );
        }
      }
    }

    inline OutputType Dirac( const InputType& iX ) const
    {
      if( vnl_math_abs( iX ) > this->m_Epsilon )
      {
        return 0.;
      }
      else
      {
        return 0.5 * this->m_InvEpsilon * ( 1.+vcl_cos( vnl_math::pi*iX * this->m_InvEpsilon ) );
      }
    }
  };
}
#endif
