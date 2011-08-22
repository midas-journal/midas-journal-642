#ifndef __itkRegionBasedLevelSetFunction_txx_
#define __itkRegionBasedLevelSetFunction_txx_

#include "itkRegionBasedLevelSetFunction.h"

namespace itk
{

template < class TInput,
  class TFeature,
  class TSharedData >
RegionBasedLevelSetFunction< TInput,
  TFeature,
  TSharedData >
::RegionBasedLevelSetFunction()
{
  m_AreaWeight = NumericTraits<ScalarValueType>::Zero;
	m_Lambda1 = 1;
	m_Lambda2 = 1;
	m_OverlapPenaltyWeight = NumericTraits<ScalarValueType>::Zero;
	m_VolumeMatchingWeight = NumericTraits<ScalarValueType>::Zero;
	m_Volume = 0;
  sharedData = 0;
  updatedC = false;
  updatedH = false;

  m_CurvatureWeight = NumericTraits<ScalarValueType>::Zero;
}

/* Computes the Heaviside function and stores it in hVals */
template < class TInput,
  class TFeature,
  class TSharedData >
void RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::ComputeHImage()
{
  // The phi function
  InputImageConstPointer contourImage = this->m_InitialImage;
  InputImagePointer hBuffer = this->sharedData->hVals[this->m_FunctionId];

  // Iterator for the phi function
  ConstImageIteratorType constIt( contourImage,
    contourImage->GetRequestedRegion() );
  ImageIteratorType It( hBuffer, hBuffer->GetRequestedRegion() );

  for( It.GoToBegin(), constIt.GoToBegin(); !constIt.IsAtEnd();
    ++It, ++constIt )
  {
    ScalarValueType hVal = m_DomainFunction->Heaviside( - constIt.Get() );
    It.Set( hVal );
  }
}

template < class TInput,
  class TFeature,
  class TSharedData >
void
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::UpdateSharedData( bool forceUpdate )
{
  if ( forceUpdate )
  {
    updatedC = false;
    updatedH = false;
  }
  if ( !updatedH )
  {
    // update H
    updatedH = true;
    this->ComputeHImage();

    // Must update all H before updating C
    return;
  }
  if ( !updatedC )
  {
    updatedC = true;
    this->ComputeParameters();
  }

  this->SpecialProcessing();

  unsigned int fId = this->m_FunctionId;

  if ( sharedData->cDens[fId] == 0 )
    sharedData->cVals[fId] = 0;
  else
    sharedData->cVals[fId] = sharedData->cNums[fId] / sharedData->cDens[fId];

	if ( sharedData->cBDen[fId] == 0 )
    sharedData->cB[fId] = 0;
  else
  	sharedData->cB[fId] = sharedData->cBNum[fId] / sharedData->cBDen[fId];

#ifndef NDEBUG
  std::cout << fId << std::endl;
	std::cout << "Background: " << std::endl;
	std::cout << this->sharedData->cBNum[fId] << ' ' << this->sharedData->cBDen[fId] << ' ';
	std::cout << this->sharedData->cB[fId] << std::endl;
	std::cout << "Foreground: " << std::endl;
	std::cout << sharedData->cNums[fId] << ' ' << this->sharedData->cDens[fId] << ' ';
	std::cout << this->sharedData->cVals[fId] << std::endl;
#endif
}

template < class TInput,
  class TFeature,
  class TSharedData >
double
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::m_WaveDT = 1.0/(2.0 * ImageDimension);

template < class TInput,
  class TFeature,
  class TSharedData >
double
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::m_DT     = 1.0/(2.0 * ImageDimension);

template < class TInput,
  class TFeature,
  class TSharedData >
typename RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >::TimeStepType
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::ComputeGlobalTimeStep(void *GlobalData) const
{
/* Computing the time-step for stable curve evolution */

  TimeStepType dt;

  ACGlobalDataStruct *d = (ACGlobalDataStruct *)GlobalData;

  if (vnl_math_abs(d->m_MaxCurvatureChange) > 0.0)
  {
    if (d->m_MaxGlobalChange > 0.0)
    {
      dt = vnl_math_min( ( this->m_WaveDT / d->m_MaxGlobalChange ),
      ( this->m_DT / d->m_MaxCurvatureChange ) );
    }
    else
    {
      dt = this->m_DT / d->m_MaxCurvatureChange;
    }
  }
  else
  {
    if (d->m_MaxGlobalChange > 0.0)
    {
//       std::cout << this->m_WaveDT << std::endl;
      dt = this->m_WaveDT / d->m_MaxGlobalChange;
    }
    else
    {
      dt = 0.0;
    }
  }

  // Reset the values
  d->m_MaxCurvatureChange   = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxGlobalChange		  = NumericTraits<ScalarValueType>::Zero;

  return dt;
}

template < class TInput,
  class TFeature, class TSharedData >
typename RegionBasedLevelSetFunction< TInput,
  TFeature, TSharedData >::
ScalarValueType
RegionBasedLevelSetFunction< TInput,
  TFeature, TSharedData >::
ComputeCurvatureTerm(
  const NeighborhoodType &itkNotUsed(neighborhood),
  const FloatOffsetType &itkNotUsed(offset), GlobalDataStruct *gd)
{
  // Calculate the mean curvature
  ScalarValueType curvature_term = NumericTraits<ScalarValueType>::Zero;
  unsigned int i, j;


  for (i = 0; i < ImageDimension; i++)
    {
    for(j = 0; j < ImageDimension; j++)
      {
      if(j != i)
        {
        curvature_term -= gd->m_dx[i] * gd->m_dx[j] * gd->m_dxy[i][j];
        curvature_term += gd->m_dxy[j][j] * gd->m_dx[i] * gd->m_dx[i];
        }
      }
    }

  return (curvature_term / gd->m_GradMagSqr );
}

template < class TInput,
  class TFeature,
  class TSharedData >
typename RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >::PixelType
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::ComputeUpdate( const NeighborhoodType &it, void *globalData,
  const FloatOffsetType& offset )
{
  // Access the neighborhood center pixel of phi
  const ScalarValueType inputValue = it.GetCenterPixel();

  ScalarValueType laplacian, laplacian_term, curvature_term, globalTerm;

  // Access the global data structure
  ACGlobalDataStruct *gd = (ACGlobalDataStruct *)globalData;

  // Compute the Hessian matrix and various other derivatives.  Some of these
  // derivatives may be used by overloaded virtual functions.
  gd->m_GradMagSqr = 1.0e-6;
  for( unsigned int i = 0 ; i < ImageDimension; i++)
  {
    const unsigned int positionA =
      static_cast< unsigned int >( this->m_Center + this->m_xStride[i] );
    const unsigned int positionB =
      static_cast< unsigned int >( this->m_Center - this->m_xStride[i] );

    gd->m_dx[i] = 0.5 * ( it.GetPixel( positionA ) - it.GetPixel( positionB ) );
    gd->m_dxy[i][i] =
      it.GetPixel( positionA ) + it.GetPixel( positionB ) - 2.0 * inputValue;
    gd->m_dx_forward[i]  = it.GetPixel( positionA ) - inputValue;
    gd->m_dx_backward[i] = inputValue - it.GetPixel( positionB );
    gd->m_GradMagSqr += gd->m_dx[i] * gd->m_dx[i];

    for( unsigned int j = i+1; j < ImageDimension; j++ )
    {
      const unsigned int positionAa = static_cast<unsigned int>(
        this->m_Center - this->m_xStride[i] - this->m_xStride[j] );
      const unsigned int positionBa = static_cast<unsigned int>(
        this->m_Center - this->m_xStride[i] + this->m_xStride[j] );
      const unsigned int positionCa = static_cast<unsigned int>(
        this->m_Center + this->m_xStride[i] - this->m_xStride[j] );
      const unsigned int positionDa = static_cast<unsigned int>(
        this->m_Center + this->m_xStride[i] + this->m_xStride[j] );

      gd->m_dxy[i][j] = gd->m_dxy[j][i] = 0.25 *(
        it.GetPixel( positionAa ) -
        it.GetPixel( positionBa ) -
        it.GetPixel( positionCa ) +
        it.GetPixel( positionDa ) );
    }
  }

  ScalarValueType dh = m_DomainFunction->Dirac( - inputValue );

  // Computing the curvature term
  // Used to regularized using the length of contour
  // What is CurvatureSpeed?
  if ( ( dh != 0. ) &&
    ( this->m_CurvatureWeight != NumericTraits< ScalarValueType >::Zero ) )
  {
    curvature_term =
      this->m_CurvatureWeight *
      this->CurvatureSpeed(it, offset) *
      this->ComputeCurvatureTerm( it, offset, gd ) *
      gd->m_GradMagSqr * dh;

    gd->m_MaxCurvatureChange =
      vnl_math_max( gd->m_MaxCurvatureChange, vnl_math_abs( curvature_term ) );
  }
  else
  {
    curvature_term = NumericTraits<ScalarValueType>::Zero;
  }

  // Computing the laplacian term
  // Used in maintaining squared distance function
  if( this->m_LaplacianSmoothingWeight != NumericTraits<ScalarValueType>::Zero)
  {
    laplacian = NumericTraits<ScalarValueType>::Zero;

    // Compute the laplacian using the existing second derivative values
    for(unsigned int i = 0;i < ImageDimension; i++)
    {
      laplacian += gd->m_dxy[i][i];
    }

	  // Use the laplacian to maintain signed distance function
    // What is LaplacianSmoothingSpeed ?
    // Why do we have 0.1 * LaplacianSmoothingWeight ?
    // Why do we have to subtract the curvature_term ?
    laplacian_term =
      0.1*( this->m_LaplacianSmoothingWeight *
      LaplacianSmoothingSpeed(it,offset, gd) * laplacian - curvature_term);
  }
  else
  {
    laplacian_term = NumericTraits<ScalarValueType>::Zero;
  }

  // Update value from curvature length and laplacian term
  PixelType updateVal =
    static_cast< PixelType >( curvature_term + laplacian_term );

  /* Compute the globalTerm - rms difference of image with c_0 or c_1*/
  if ( dh != 0. )
    globalTerm = dh * this->computeGlobalTerm( inputValue, it.GetIndex() );
  else
    globalTerm = 0;

  /* Final update value is the local terms of curvature lengths and laplacian
  squared distances - global terms of rms differences of image and piecewise
  constant regions*/
  updateVal = updateVal - globalTerm;

  /* If MaxGlobalChange recorded is lower than the current globalTerm */
  if(vnl_math_abs( gd->m_MaxGlobalChange) < vnl_math_abs( globalTerm ) )
  {
    gd->m_MaxGlobalChange = globalTerm;
  }

  return updateVal;
}

/* Computes the fidelity term (eg: (intensity - mean)2 ).
Most of the code is concerned with using the appropriate combination
of Heaviside and dirac delta for each part of the fidelity term.
- the final dH is the dirac delta term corresponding to the current
level set we are updating. */
template < class TInput, class TFeature, class TSharedData >
typename RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::ScalarValueType
RegionBasedLevelSetFunction< TInput, TFeature, TSharedData >
::computeGlobalTerm(
const ScalarValueType& inputPixel,
const InputIndexType& inputIndex )
{
  unsigned int fId = this->m_FunctionId;
  unsigned int pr = 1; // computes if it belongs to background
  unsigned int s = 0; // accumulates the overlap across all functions

  // Assuming only 1 level set function to be present
  FeatureIndexType featIndex = static_cast< FeatureIndexType >( inputIndex );

  const FeaturePixelType featureVal =
    this->m_FeatureImage->GetPixel ( inputIndex );

  ScalarValueType globalTerm = 0;

  ScalarValueType overlapTerm = 0.;
  // This conditional statement computes the amount of overlap s
  // and the presence of background pr
  if ( this->sharedData->FunctionCount > 1 )
  {
    featIndex = this->sharedData->GetFeatureIndex( fId, inputIndex );
    computeOverlapParameters( featIndex, s, pr );
  }

  ScalarValueType inTerm = this->computeInternalTerm( featureVal, featIndex, fId );
  ScalarValueType outTerm = this->m_Lambda2 * this->computeExternalTerm( featureVal, featIndex, pr );
  overlapTerm = this->computeOverlapTerm( s );
  ScalarValueType regularizationTerm = 2 * this->m_VolumeMatchingWeight *
    ( sharedData->cDens[fId] - this->m_Volume );
  //regularizationTerm -= this->m_Nu;
  //NOTE: regularizationTerm here MUST take into account the curvature term!!!

  globalTerm = - inTerm + outTerm - overlapTerm - regularizationTerm;

  return globalTerm;
}

} // end namespace

#endif
