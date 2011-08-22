#ifndef __itkDenseMultiphaseLevelSetImageFilter_h_
#define __itkDenseMultiphaseLevelSetImageFilter_h_

#include "itkMultiphaseDenseSegmentationFiniteDifferenceImageFilter.h"
#include "itkScalarChanAndVeseSharedFunctionData.h"
#include "itkRegionOfInterestImageFilter.h"

#ifdef VIZU
  #include "vtkVisualize2DImplicitFunction.h"
  #include "vtkVisualize3DImplicitFunction.h"
#endif

namespace itk
{
template < class TInput,
  class TFeature,
  class TFunction,
  typename TOutputPixel = typename TInput::PointType::CoordRepType,
  class TSharedData =
    ScalarChanAndVeseSharedFunctionData<
    TInput, TFeature > >
class ITK_EXPORT DenseMultiphaseLevelSetImageFilter:
public DenseSegmentationFiniteDifferenceImageFilter<
  TInput,
  Image< TOutputPixel, TInput::ImageDimension >,
  TFunction >
{
public:

  /** Repeat definition from Superclass to satisfy Borland compiler quirks */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInput::ImageDimension );

  typedef Image< TOutputPixel, ImageDimension > OutputImageType;

  typedef DenseMultiphaseLevelSetImageFilter Self;
  typedef DenseSegmentationFiniteDifferenceImageFilter<TInput,
    OutputImageType, TFunction > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DenseMultiphaseLevelSetImageFilter,
    MultiphaseDenseSegmentationFiniteDifferenceImageFilter );

  void PrintSelf( std::ostream& os, Indent indent ) const;

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
    vtkVisualize2DImplicitFunctionType;
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
    this->ProcessObject::SetNthInput( 0, const_cast< FeatureImageType * >(f) );
  }

  virtual FeatureImageType * GetFeatureImage()
  {
    return (static_cast< FeatureImageType*>(this->ProcessObject::GetInput(0)));
  }

protected:
  DenseMultiphaseLevelSetImageFilter()
    { sharedData = SharedDataType::New(); }
  ~DenseMultiphaseLevelSetImageFilter(){}

#if VIZU
  vtkVisualize2DImplicitFunctionType func;
#endif

  SharedDataPointer sharedData;

  virtual void Initialize();
  virtual void InitializeIteration();
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDenseMultiphaseLevelSetImageFilter.txx"
#endif

#endif
