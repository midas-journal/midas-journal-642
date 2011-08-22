#ifndef __itkMultiLevelSetTracking_txx
#define __itkMultiLevelSetTracking_txx

#include "itkMultiLevelSetTracking.h"

namespace itk
{
template < class TInputImage,
  class TFeatureImage,
  class TInternalImage,
  class TLevelSetFunction,
  class TMultiPhaseFilter >
MultiLevelSetTracking< TInputImage, TFeatureImage, TInternalImage,
TLevelSetFunction, TMultiPhaseFilter >::
MultiLevelSetTracking()
{
  m_ROIfilter = ROIFilterType::New();
  m_LevelSetFilter = MultiLevelSetType::New();
  m_CastOutput = OutputCastType::New();

  m_Iterations = 10;
  m_RMSError = 0;
  m_Epsilon = 1;
  m_CurvatureWeight = 0;
  m_AreaWeight = 0;
  m_Lambda1 = 1;
  m_Lambda2 = 1;
  m_OverlapPenaltyWeight = 0;
  m_LaplacianSmoothingWeight = 0;
  m_VolumeMatchingWeight = 0;
  m_Volume = 100;
  m_FunctionCount = 1;
  m_LargestCellRadius = 8.0;

  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->Superclass::SetNumberOfRequiredOutputs( 1 );

  this->Superclass::SetNthOutput( 0, TInputImage::New() );
}

template < class TInputImage,
  class TFeatureImage,
  class TInternalImage,
  class TLevelSetFunction,
  class TMultiPhaseFilter >
typename MultiLevelSetTracking< TInputImage,TFeatureImage,
TInternalImage, TLevelSetFunction, TMultiPhaseFilter >::ImagePointer
MultiLevelSetTracking< TInputImage,TFeatureImage,
TInternalImage, TLevelSetFunction, TMultiPhaseFilter >::
ComputeLevelSetFunction( ImagePointer input, unsigned int index )
{
  // Allocate an image of same size as input and then assign it to the
  typename ThresholdFilterType::Pointer thresh = ThresholdFilterType::New();;
  thresh->SetLowerThreshold( index );
  thresh->SetUpperThreshold( index );
  thresh->SetInsideValue( 1 );
  thresh->SetOutsideValue( 0 );
  thresh->SetInput( input );
  thresh->Update();

  typename MaurerType::Pointer maurer = MaurerType::New();
  maurer->SetInput( thresh->GetOutput() );
  maurer->SetSquaredDistance( 0 );
  maurer->SetUseImageSpacing( 1 );
  maurer->SetInsideIsPositive( 0 );
  maurer->Update();

  return ( maurer->GetOutput() );
}

template < class TInputImage,
  class TFeatureImage,
  class TInternalImage,
  class TLevelSetFunction,
  class TMultiPhaseFilter >
typename MultiLevelSetTracking< TInputImage,TFeatureImage,
TInternalImage, TLevelSetFunction, TMultiPhaseFilter >::
ImagePointer
MultiLevelSetTracking< TInputImage,TFeatureImage,
TInternalImage, TLevelSetFunction, TMultiPhaseFilter >::
Resample
( ImagePointer iInput, ImageSpacingType iSpacing, ImageSizeType iSize,
ImagePointType iOrigin )
{
  // create the resample filter, transform and interpolator
  TransformPointer transform = TransformType::New();
  transform->SetIdentity();

  InterpolatorPointer interp = InterpolatorType::New();

  ResampleFilterPointer m_resamp = ResampleFilterType::New();
  m_resamp->SetTransform( transform );
  m_resamp->SetInterpolator( interp );
  m_resamp->SetInput( iInput );
  m_resamp->SetSize( iSize );
  m_resamp->SetOutputOrigin( iOrigin );
  m_resamp->SetOutputSpacing( iSpacing );
  m_resamp->SetDefaultPixelValue( 0 );
  m_resamp->Update();

  return(m_resamp->GetOutput());
}

template < class TInputImage,
  class TFeatureImage,
  class TInternalImage,
  class TLevelSetFunction,
  class TMultiPhaseFilter >
void
MultiLevelSetTracking< TInputImage,TFeatureImage, TInternalImage,
TLevelSetFunction, TMultiPhaseFilter >::
GenerateData()
{
  HeavisideFunctionType Heaviside( m_Epsilon );

  ImagePointer input;
  {
    typename InputCastType::Pointer m_CastInput =
      InputCastType::New();
    m_CastInput->SetInput( this->GetInput() );
    m_CastInput->Update();
    input = m_CastInput->GetOutput();
    input->DisconnectPipeline();
  }

  ImageSpacingType spacing = input->GetSpacing();
  ImagePointType origin = input->GetOrigin();
  ImageSizeType inputSize = input->GetLargestPossibleRegion().GetSize();

  ImageSpacingType  subSpacing;
  ImageSizeType     subSize;
  for(unsigned int i = 0; i < ImageDimension; i++)
  {
    subSpacing[i] = spacing[i] * m_Sampling[i];
    subSize[i]    = static_cast<ImageSizeValueType>(inputSize[i]/m_Sampling[i]);
  }

  input = Resample( input, subSpacing, subSize, origin );

  // Find the number of cells and their lookup table
  VectorType m_Lookup;
  typename SampleType::Pointer centroidSample = SampleType::New();
  centroidSample->SetMeasurementVectorSize( ImageDimension );
  {
    // convert the image into a collection of objects
    typename ShapeConverterType::Pointer shapeConverter =
      ShapeConverterType::New();
    shapeConverter->SetInput( input );
    shapeConverter->SetBackgroundValue( 0 );
    shapeConverter->Update();
    ShapeLabelMapPointer shapeLabelMap = shapeConverter->GetOutput();

    // Get the number of label objects
    m_FunctionCount = shapeLabelMap->GetNumberOfLabelObjects();

    // Computing the 0-indexed lookup table
    LabelObjectContainerType container =
      shapeLabelMap->GetLabelObjectContainer();

    unsigned int label = 0;
    m_Lookup.resize( m_FunctionCount );
    for(  LabelObjectIterator it = container.begin();
      it != container.end(); it++, label++ )
    {
      m_Lookup[label] = static_cast< unsigned int >( it->first );
    }

    CentroidVectorType centroid;
    InputPointType idx;
    for( unsigned short i = 0; i < m_FunctionCount; i++)
    {
      idx = shapeLabelMap->GetLabelObject( m_Lookup[i] )->GetCentroid();
      for( unsigned int j = 0; j < ImageDimension; j++)
        centroid[j] = static_cast<float>( idx[j] );
      centroidSample->PushBack( centroid );
    }
  }
#ifndef NDEBUG
  std::cout << "Total cells being tracked: " << m_FunctionCount << std::endl;
#endif
  m_LevelSetFilter->SetFunctionCount( m_FunctionCount );
  m_LevelSetFilter->SetLookup( m_Lookup );

  ImagePointer featureImage;
  {
    typename FeatureCastType::Pointer m_CastFeature =
      FeatureCastType::New();
    m_CastFeature->SetInput( m_Feature );
    m_CastFeature->Update();
    featureImage = Resample( m_CastFeature->GetOutput(), subSpacing,
      subSize, origin );
    featureImage->DisconnectPipeline();
  }

  m_LevelSetFilter->SetFeatureImage( featureImage );
  m_Feature = 0; // Free the memory by de-referencing the smart pointer

  // Compute the bounding box and use ROI filter
  m_ROIfilter->SetInput( input );

  InputIndexType cellExtent;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    cellExtent[j] =
      static_cast<InputIndexValueType>( m_LargestCellRadius/subSpacing[j] );

  for( unsigned short i = 0; i < m_FunctionCount; i++)
  {
    // Code to calculate ROIs
    InputIndexType start,end;
    InputSizeType size;

    // This is in physical coordinates
    InputIndexValueType idx;

    CentroidVectorType centroid = centroidSample->GetMeasurementVector(i);
    for( unsigned int j = 0; j < ImageDimension; j++ )
    {
      idx = static_cast< InputIndexValueType >( centroid[j]/subSpacing[j] );

      if ( idx >= cellExtent[j] )
        start[j] = static_cast< InputIndexValueType >( idx - cellExtent[j]);
      else
        start[j] = 0;

      if ( static_cast< InputSizeValueType >( idx + cellExtent[j] )
        < subSize[j] )
        end[j] = static_cast< InputIndexValueType>( idx + cellExtent[j] );
      else
        end[j] = static_cast< InputIndexValueType>( subSize[j]-1 );

      size[j] = end[j] - start[j];
    }

    InputRegionType region;
    region.SetSize( size );
    region.SetIndex( start );

    m_ROIfilter->SetRegionOfInterest( region );
    m_ROIfilter->Update();
    // Note that the origin of the roi image is in physical coordinates
    // GetOutput() returns an image strictly defined only inside the region

    // Compute their distance maps in a loop only in the requested region
    ImagePointer contourImage = ComputeLevelSetFunction(
      m_ROIfilter->GetOutput(), m_Lookup[i] );

#ifndef NDEBUG
    std::cout << "Label: " << i << " Origin: " << contourImage->GetOrigin() << std::endl;
    std::cout << "Size: " << contourImage->GetLargestPossibleRegion().GetSize() << std::endl;
#endif

    m_LevelSetFilter->SetLevelSet( i, contourImage );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetDomainFunction( &Heaviside );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetCurvatureWeight( m_CurvatureWeight );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetAreaWeight( m_AreaWeight );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetLambda1( m_Lambda1 );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetLambda2( m_Lambda2 );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetOverlapPenaltyWeight( m_OverlapPenaltyWeight );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetLaplacianSmoothingWeight( m_LaplacianSmoothingWeight );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetVolumeMatchingWeight( m_VolumeMatchingWeight );
    m_LevelSetFilter->GetDifferenceFunction(i)->SetVolume( m_Volume );
  }
  m_ROIfilter = 0;

  if ( m_FunctionCount > 20 )
  {
    TreePointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample( centroidSample );
    treeGenerator->SetBucketSize( 5 );
    treeGenerator->Update();
    KdTreePointer kdTree = treeGenerator->GetOutput();

    m_LevelSetFilter->SetKdTree( kdTree );
  }

  // Initialize the levelset filter
  m_LevelSetFilter->SetNumberOfIterations( m_Iterations );
  m_LevelSetFilter->SetMaximumRMSError( m_RMSError );
  m_LevelSetFilter->SetUseImageSpacing( 1 );
  m_LevelSetFilter->Update();

  featureImage = 0; // Delete the featureImage
  ImagePointer output = Resample( m_LevelSetFilter->GetOutput(), spacing,
    inputSize, origin );

  // Allocate an output image
  m_CastOutput->SetInput( output );
  m_CastOutput->GraftOutput( this->GetOutput() );
  m_CastOutput->Update();
  std::cout << "Finished processing..." << std::endl;

  this->GraftOutput( m_CastOutput->GetOutput() );
}

template < class TInputImage,
  class TFeatureImage,
  class TInternalImage,
  class TLevelSetFunction,
  class TMultiPhaseFilter >
void
MultiLevelSetTracking< TInputImage,TFeatureImage, TInternalImage,
TLevelSetFunction, TMultiPhaseFilter >::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Name Of Class: " << GetNameOfClass() << std::endl;
  os << indent << "m_Iterations: " << GetIterations() << std::endl;
  os << indent << "m_RMSError: " << GetRMSError() << std::endl;
  os << indent << "m_Epsilon: " << this->m_Epsilon << std::endl;
  os << indent << "m_CurvatureWeight: " << this->m_CurvatureWeight << std::endl;
  os << indent << "m_AreaWeight: "<< this->m_AreaWeight << std::endl;
  os << indent << "m_Lambda1: " << this->m_Lambda1 << std::endl;
  os << indent << "m_Lambda2: " << this->m_Lambda2 << std::endl;
  os << indent << "m_OverlapPenaltyWeight: " << this->m_OverlapPenaltyWeight << std::endl;
  os << indent << "m_LaplacianSmoothingWeight: " << this->m_LaplacianSmoothingWeight << std::endl;
  os << indent << "m_VolumeMatchingWeight: " << this->m_VolumeMatchingWeight << std::endl;
  os << indent << "m_Volume: " << this->m_Volume << std::endl;
  os << indent << "m_FunctionCount: " << GetFunctionCount()<< std::endl;
  os << indent << "m_LargestCellRadius: " << GetLargestCellRadius()<< std::endl;
}

} /* end namespace itk */

#endif
