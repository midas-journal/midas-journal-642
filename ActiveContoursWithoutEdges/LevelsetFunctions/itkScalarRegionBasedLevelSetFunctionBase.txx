#ifndef __itkScalarRegionBasedLevelSetFunctionBase_txx_
#define __itkScalarRegionBasedLevelSetFunctionBase_txx_

#include "itkScalarRegionBasedLevelSetFunctionBase.h"

namespace itk {

/* Calculates the numerator and denominator for c_i for each region. As part of
the optimization, it is called once at the beginning of the code, and then the
cNum and cDen are updated during the evolution without iterating through the
entire image. */
template < class TInputImage, class TFeatureImage, class TSharedData >
void
ScalarRegionBasedLevelSetFunctionBase< TInputImage, TFeatureImage, TSharedData >
::ComputeParameters()
{
  unsigned int fId = this->m_FunctionId;

  this->sharedData->cDens[fId] = 0;
  this->sharedData->cNums[fId] = 0;
  this->sharedData->cVals[fId] = 0;
  this->sharedData->cBDen[fId] = 0;
  this->sharedData->cBNum[fId] = 0;
  this->sharedData->cB[fId] = 0;

	FeatureImageConstPointer featureImage = this->m_FeatureImage;

  ImageIteratorType It( this->sharedData->hVals[fId],
    this->sharedData->hVals[fId]->GetLargestPossibleRegion() );
  ConstFeatureIteratorType fIt( this->m_FeatureImage,
    this->m_FeatureImage->GetLargestPossibleRegion() );

  FeaturePixelType featureVal;
  FeatureIndexType globalIndex;
  InputIndexType itInputIndex, inputIndex;
  InputPixelType hVal, inputVal;
  ListPixelType L;

  for( It.GoToBegin(), fIt.GoToBegin(); !It.IsAtEnd();
    ++It, ++fIt )
  {
    inputVal = It.Get();
    featureVal = fIt.Get();
    inputIndex = It.GetIndex();

		globalIndex = this->sharedData->GetFeatureIndex( fId, inputIndex );

    L = this->sharedData->lImage->GetPixel( globalIndex );

    bool inBgrnd = true; // assume the pixel is in background
    for( ListPixelConstIterator it = L.begin(); it != L.end(); ++it )
    {
      itInputIndex = this->sharedData->GetIndex( *it, globalIndex );
      hVal = this->sharedData->hVals[*it]->GetPixel( itInputIndex );

      if ( hVal > 0.5 )
      {
        // inside the level-set function
        inBgrnd = false;

        if (*it == fId)
        {
            this->sharedData->cNums[fId] += featureVal;
            this->sharedData->cDens[fId] += 1.;
        }
      }
    }

    // if the pixel belongs to the background
    if ( inBgrnd )
    {
      this->sharedData->cBNum[fId] += featureVal;
      this->sharedData->cBDen[fId] += 1.;
    }
  }
}


template < class TInputImage, class TFeatureImage, class TSharedData >
void
ScalarRegionBasedLevelSetFunctionBase< TInputImage, TFeatureImage, TSharedData >
::computeOverlapParameters( const FeatureIndexType globalIndex, unsigned int& s, unsigned int& pr )
{
// This conditional statement computes the amount of overlap s
// and the presence of background pr
  unsigned int fId = this->m_FunctionId;

  ListPixelType L;
  L = this->sharedData->lImage->GetPixel( globalIndex );

  InputPixelType hVal;
  InputIndexType otherIndex;

  for ( ListPixelIterator it = L.begin(); it != L.end(); ++it )
  {
    unsigned int id = *it;
    if ( id != fId )
    {
      otherIndex = this->sharedData->GetIndex( id, globalIndex );
      hVal = this->sharedData->hVals[id]->GetPixel( otherIndex );
      if ( hVal > 0.5 )
      {
        s ++;
        pr = 0;
      }
    }
  }
}

/* Performs the narrow-band update of the Heaviside function for each voxel. The
characteristic function of each region is recomputed (note the shared
data which contains information from the other level sets). Using the
new H values, the previous c_i are updated. */
template < class TInputImage, class TFeatureImage, class TSharedData >
void
ScalarRegionBasedLevelSetFunctionBase< TInputImage, TFeatureImage, TSharedData >
::UpdatePixel( const unsigned int& idx, NeighborhoodIterator< TInputImage >
&iterator, InputPixelType &newValue, bool &status )
{
  unsigned int fId = this->m_FunctionId;

  // For each affected h val: h val = new hval (this will dirty some cvals)
  InputIndexType inputIndex = iterator.GetIndex( idx );
  FeatureIndexType globalIndex =
    this->sharedData->GetFeatureIndex( fId, inputIndex );

  FeaturePixelType featureVal =
    this->m_FeatureImage->GetPixel( inputIndex );

  ScalarValueType oldH = this->sharedData->hVals[fId]->GetPixel( inputIndex );
  ScalarValueType newH = this->m_DomainFunction->Heaviside( - newValue );

  // Check if it is in other foreground
  ListPixelType L = this->sharedData->lImage->GetPixel( globalIndex );
  InputIndexType itInputIndex;
  ScalarValueType hVal;

  bool inBgrnd = true; // assume the pixel is in background
  for( ListPixelType::const_iterator it = L.begin(); it != L.end(); ++it )
  {
    itInputIndex = this->sharedData->GetIndex( *it, globalIndex );
    hVal = this->sharedData->hVals[*it]->GetPixel( itInputIndex );

    if ( ( hVal > 0.5 ) && ( *it != fId ) )
    {
      inBgrnd = false; // belongs to foreground elsewhere
    }
  }

  // if pixel belonged to current foreground but not anymore so
  if ( ( oldH > 0.5 ) && ( newH <= 0.5 ) )
  {
    this->sharedData->cDens[fId]--;
    this->sharedData->cNums[fId] -= featureVal;

    // have to update level-set backgrounds overlapping
    // at the current pixel
    if ( inBgrnd )
    {
      for( ListPixelType::const_iterator it = L.begin(); it != L.end(); ++it )
      {
        this->sharedData->cBNum[*it] += featureVal;
        this->sharedData->cBDen[*it]++;
      }
    }
  }

  // if pixel entered the foreground
  if ( ( oldH <= 0.5 ) & ( newH > 0.5 ) )
  {
    this->sharedData->cDens[fId]++;
    this->sharedData->cNums[fId] += featureVal;

    // have to update level-set backgrounds overlapping
    // at the current pixel
    if ( inBgrnd )
    {
      for( ListPixelType::const_iterator it = L.begin(); it != L.end(); ++it )
      {
        this->sharedData->cBNum[*it] -= featureVal;
        this->sharedData->cBDen[*it]--;
      }
    }
  }

  this->sharedData->hVals[fId]->SetPixel( inputIndex, newH );
}


} // end namespace itk
#endif
