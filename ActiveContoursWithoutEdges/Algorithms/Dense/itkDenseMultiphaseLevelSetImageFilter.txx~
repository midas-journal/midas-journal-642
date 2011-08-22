#ifndef __itkDenseMultiphaseLevelSetImageFilter_txx
#define __itkDenseMultiphaseLevelSetImageFilter_txx

#include "itkDenseMultiphaseLevelSetImageFilter.h"
#include "itkScalarChanAndVeseSharedFunctionData.h"

namespace itk
{
template < class TInput,
  class TFeature,
  class TFunction,
  class TOutputPixel,
  class TSharedData >
void
DenseMultiphaseLevelSetImageFilter< TInput, TFeature,
  TFunction, TOutputPixel, TSharedData >::
Initialize()
{
  // Set the feature image for the individual level-set functions
  for( unsigned int i = 0; i < this->m_FunctionCount; i++)
  {
    InputImagePointer input = this->m_LevelSet[i];
    InputPointType origin = input->GetOrigin();
    InputSpacingType spacing = input->GetSpacing();

    // In the context of the global coordinates
    FeatureIndexType start;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      start[j] = static_cast<FeatureIndexValueType>( origin[j]/spacing[j] );

    // Defining roi region
    FeatureRegionType region;
    region.SetSize( input->GetLargestPossibleRegion().GetSize() );
    region.SetIndex( start );

    // Initialize the ROI filter with the feature image
    ROIFilterPointer roi = ROIFilterType::New();
    roi->SetInput( this->GetFeatureImage() );
    roi->SetRegionOfInterest( region );
    roi->Update();

    // Assign roi output
    FeatureImagePtr feature = roi->GetOutput();
    this->m_DifferenceFunctions[i]->SetFeatureImage( feature );
    this->m_DifferenceFunctions[i]->SetInitialImage( input );
  }

#ifndef NDEBUG
  std::cout << "Set feature image regions of interest " << std::endl;
#endif

  // Initialize the function count in sharedData
  sharedData->SetFunctionCount ( this->m_FunctionCount );

  // Set the KdTree pointer
  if ( this->m_KdTree )
  {
    sharedData->SetKdTree( this->m_KdTree );
  }

  for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
  {
    this->m_DifferenceFunctions[i]->SetFunctionId( i );

    sharedData->CreateHVals (
      i, this->m_LevelSet[i]->GetSpacing(),
      this->m_LevelSet[i]->GetOrigin(),
      this->m_LevelSet[i]->GetLargestPossibleRegion() );

    // Share the sharedData structure
    this->m_DifferenceFunctions[i]->SetSharedData( sharedData );
  }

  sharedData->AllocateListImage(
    this->GetFeatureImage()->GetLargestPossibleRegion(),
    this->GetFeatureImage()->GetSpacing() );

#ifndef NDEBUG
  std::cout << "Allocated list image " << std::endl;
#endif

  sharedData->PopulateListImage();

#ifndef NDEBUG
  std::cout << "Populated list image " << std::endl;
#endif

  Superclass::Initialize();

  for (unsigned int i = 0; i < this->m_FunctionCount; i++)
  {
    this->m_DifferenceFunctions[i]->UpdateSharedData(true);
  }

  for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
  {
    this->m_DifferenceFunctions[i]->UpdateSharedData( false );
  }

#ifdef VIZU
  func.SetFunctionCount( this->m_FunctionCount );
  func.SetFeatureImage( this->GetFeatureImage() );
#endif
}

/** Overrides parent implementation */
// This function is called at the end of each iteration
template < class TInput,
  class TFeature,
  class TFunction,
  class TOutputPixel,
  class TSharedData >
void
DenseMultiphaseLevelSetImageFilter< TInput, TFeature,
  TFunction, TOutputPixel, TSharedData > ::
InitializeIteration()
{
  Superclass::InitializeIteration();

  for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
  {
    this->m_DifferenceFunctions[i]->SetInitialImage( this->m_LevelSet[i] );
    this->m_DifferenceFunctions[i]->UpdateSharedData ( true );
  }

  for ( unsigned int i = 0; i < this->m_FunctionCount; i++ )
  {
    this->m_DifferenceFunctions[i]->UpdateSharedData ( false );

    //Visualization routine
#ifdef VIZU
    func.SetInputImage( this->m_LevelSet[i] );
    func.Visualize( this->m_ElapsedIterations, 0, i );
#endif
  }

  // Estimate the progress of the filter
  this->SetProgress ( ( ( float ) this->m_ElapsedIterations
    / this->m_NumberOfIterations ) );
}

template < class TInput,
  class TFeature,
  class TFunction,
  class TOutputPixel,
  class TSharedData >
void
DenseMultiphaseLevelSetImageFilter< TInput, TFeature,
  TFunction, TOutputPixel, TSharedData > ::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Class Name: " << GetNameOfClass( ) << std::endl;
}


} /* end namespace itk */

#endif
