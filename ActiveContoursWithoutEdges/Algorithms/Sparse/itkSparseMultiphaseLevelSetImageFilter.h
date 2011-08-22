#ifndef __itkSparseMultiphaseLevelSetImageFilter_h_
#define __itkSparseMultiphaseLevelSetImageFilter_h_

#include "itkMultiphaseSparseFieldLevelSetImageFilter.h"
#include "itkScalarChanAndVeseSharedFunctionData.h"
#include "itkRegionOfInterestImageFilter.h"

#ifdef VIZU
  #include "vtkVisualize2DImplicitFunction.h"
  #include "vtkVisualize3DImplicitFunction.h"
#endif

namespace itk
{
template < class TInput, // LevelSetImageType
  class TFeature,        // FeatureImageType
  class TFunction,
  typename TOutputPixel = typename TInput::PointType::CoordRepType,
  class TSharedData =
    ScalarChanAndVeseSharedFunctionData<
    TInput, TFeature > >
class ITK_EXPORT SparseMultiphaseLevelSetImageFilter:
public MultiphaseSparseFieldLevelSetImageFilter<
  TInput,
  Image< TOutputPixel, TInput::ImageDimension >,
  TFunction >
{
public:

  /** Repeat definition from Superclass to satisfy Borland compiler quirks */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInput::ImageDimension );

  typedef Image< TOutputPixel, ImageDimension > OutputImageType;

  typedef SparseMultiphaseLevelSetImageFilter Self;
  typedef MultiphaseSparseFieldLevelSetImageFilter< TInput,
    OutputImageType, TFunction > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SparseMultiphaseLevelSetImageFilter,
    MultiphaseSparseFieldLevelSetImageFilter );

  void PrintSelf( std::ostream& os, Indent indent) const;

  /** Inherited typedef from the superclass. */
  typedef TFeature FeatureImageType;
  typedef typename FeatureImageType::Pointer FeatureImagePtr;
  typedef typename FeatureImageType::PixelType FeaturePixelType;
  typedef typename FeatureImageType::IndexType FeatureIndexType;
  typedef typename FeatureIndexType::IndexValueType FeatureIndexValueType;
  typedef typename FeatureImageType::RegionType FeatureRegionType;

  /** Output image type typedefs */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::InputImagePointer InputImagePointer;
  typedef typename Superclass::InputPointType InputPointType;
  typedef typename Superclass::InputSpacingType InputSpacingType;

  typedef typename OutputImageType::ValueType ValueType;
  typedef typename OutputImageType::IndexType IndexType;

  typedef typename Superclass::TimeStepType   TimeStepType;
  typedef typename Superclass::FiniteDifferenceFunctionType
    FiniteDifferenceFunctionType;

  typedef TFunction FunctionType;
  typedef typename FunctionType::Pointer FunctionPtr;

  typedef TSharedData SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

  typedef RegionOfInterestImageFilter< FeatureImageType,
    FeatureImageType > ROIFilterType;
  typedef typename ROIFilterType::Pointer ROIFilterPointer;

#ifdef VIZU
  typedef vtkVisualize2DImplicitFunction< InputImageType, FeatureImageType >
    vtkVisualizeImplicitFunctionType;
#endif

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<TOutputPixel>) );
  /** End concept checking */
#endif

  void SetFunctionCount( unsigned int n )
  {
    Superclass::SetFunctionCount( n );
  }

  /** Set/Get the feature image to be used for speed function of the level set
   *  equation.  Equivalent to calling Set/GetInput(1, ..) */
  virtual void SetFeatureImage(const FeatureImageType *f)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< FeatureImageType *>(f) );
  }

  virtual FeatureImageType * GetFeatureImage()
  {
    return (static_cast< FeatureImageType*>(this->ProcessObject::GetInput(0)));
  }

protected:
  SparseMultiphaseLevelSetImageFilter()
  {
    this->SetNumberOfLayers(5); // Narrow-band usage
    sharedData = SharedDataType::New();
  }

  ~SparseMultiphaseLevelSetImageFilter(){}

#ifdef VIZU
  vtkVisualizeImplicitFunctionType func;
#endif
  SharedDataPointer sharedData;

  virtual void Initialize();
  virtual void InitializeIteration();
  virtual void UpdatePixel( unsigned int functionIndex,
    unsigned int idx, NeighborhoodIterator< OutputImageType > &iterator,
    ValueType &newValue, bool &status );
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSparseMultiphaseLevelSetImageFilter.txx"
#endif

#endif
