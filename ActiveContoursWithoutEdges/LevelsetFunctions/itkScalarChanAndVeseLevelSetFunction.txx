#ifndef __itkScalarChanAndVeseLevelSetFunction_txx_
#define __itkScalarChanAndVeseLevelSetFunction_txx_

#include "itkScalarChanAndVeseLevelSetFunction.h"

namespace itk {

/* Calculates the numerator and denominator for c_i for each region. As part of
the optimization, it is called once at the beginning of the code, and then the
cNum and cDen are updated during the evolution without iterating through the
entire image. */
template < class TInputImage, class TFeatureImage, class TSharedData >
void
ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ComputeParameters()
{
  Superclass::ComputeParameters();
}


template < class TInputImage, class TFeatureImage, class TSharedData >
typename ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ScalarValueType
ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
computeInternalTerm( const FeaturePixelType& iValue,
  const FeatureIndexType& iIdx, const unsigned int& fId )
{
  ScalarValueType cVals = this->sharedData->cVals[fId];
  ScalarValueType t = ( iValue - cVals );
  return this->m_Lambda1 * t * t;
}


template < class TInputImage, class TFeatureImage, class TSharedData >
typename ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ScalarValueType
ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::computeExternalTerm( const FeaturePixelType& iValue,
    const FeatureIndexType& iIdx, const unsigned int& pr )
{
  unsigned int fId = this->m_FunctionId;
  double cBgrnd = this->sharedData->cB[fId]; // bgrnd

  return pr * ( iValue - cBgrnd ) * ( iValue - cBgrnd );
}


/* Performs the narrow-band update of the Heaviside function for each voxel. The
characteristic function of each region is recomputed (note the shared
data which contains information from the other level sets). Using the
new H values, the previous c_i are updated. */
template < class TInputImage, class TFeatureImage, class TSharedData >
void
ScalarChanAndVeseLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::UpdatePixel( const unsigned int& idx, NeighborhoodIterator< TInputImage >
&iterator, ScalarValueType &newValue, bool &status )
{
  Superclass::UpdatePixel(idx,iterator,newValue,status);
}


} // end namespace itk
#endif
