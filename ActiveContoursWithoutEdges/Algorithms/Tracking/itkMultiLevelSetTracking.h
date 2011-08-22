#ifndef __itkMultiLevelSetTracking_h_
#define __itkMultiLevelSetTracking_h_

#include "itkScalarChanAndVeseLevelSetFunction.h"
#include "itkDenseMultiphaseLevelSetImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkLabelObject.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include <vnl/vnl_vector.h>

#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"

#ifdef VIZU
  #include "vtkVisualize2DImplicitFunction.h"
  #include "vtkVisualize3DImplicitFunction.h"
#endif

namespace itk {

template < class TInputImage, // LevelSetImageType
  class TFeatureImage,        // FeatureImageType
  class TInternalImage = Image< float, TInputImage::ImageDimension >,
  class TLevelSetFunction =
    ScalarChanAndVeseLevelSetFunction< TInternalImage,
    TInternalImage >,
  class TMultiPhaseFilter =
    DenseMultiphaseLevelSetImageFilter< TInternalImage,
    TInternalImage, TLevelSetFunction, float > >
class ITK_EXPORT MultiLevelSetTracking : public ImageToImageFilter<
  TInputImage, TInputImage >
{
public:

  typedef MultiLevelSetTracking Self;
  typedef ImageToImageFilter< TInputImage,TInputImage > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiLevelSetTracking, ImageToImageFilter );

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

//   typedef Image< float,ImageDimension >    ImageType;
  typedef TInternalImage                   ImageType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::PixelType    ImagePixelType;
  typedef typename ImageType::RegionType   ImageRegionType;
  typedef typename ImageType::SizeType     ImageSizeType;
  typedef typename ImageSizeType::SizeValueType ImageSizeValueType;
  typedef typename ImageType::SpacingType  ImageSpacingType;
  typedef typename ImageType::IndexType    ImageIndexType;
  typedef typename ImageIndexType::IndexValueType ImageIndexValueType;
  typedef typename ImageType::PointType    ImagePointType;

  /** Input image typedefs */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename InputImageType::IndexType InputIndexType;
  typedef typename InputIndexType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename InputSizeType::SizeValueType InputSizeValueType;
  typedef typename InputImageType::SpacingType InputSpacingType;
  typedef typename InputImageType::PointType InputPointType;
  typedef typename InputImageType::RegionType InputRegionType;

  /** Feature image typedefs */
  typedef TFeatureImage FeatureImageType;
  typedef typename FeatureImageType::Pointer FeatureImagePointer;

  /** Output image typedefs */
  typedef TInputImage OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  typedef CastImageFilter< FeatureImageType, ImageType > FeatureCastType;
  typedef CastImageFilter< InputImageType, ImageType > InputCastType;

  typedef ImageRegionConstIterator< ImageType > IteratorType;
  typedef BinaryThresholdImageFilter< ImageType, ImageType >
    ThresholdFilterType;
  typedef SignedMaurerDistanceMapImageFilter< ImageType, ImageType >
    MaurerType;
  typedef RegionOfInterestImageFilter< ImageType, ImageType >
    ROIFilterType;

  typedef ResampleImageFilter<ImageType,ImageType > ResampleFilterType;
  typedef typename ResampleFilterType::Pointer ResampleFilterPointer;
  typedef NearestNeighborInterpolateImageFunction< ImageType > InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;
  typedef AffineTransform< double, ImageDimension > TransformType;
  typedef typename TransformType::Pointer TransformPointer;

  typedef TLevelSetFunction LevelSetFunctionType;
  typedef typename LevelSetFunctionType::Pointer LevelSetFunctionPointer;

  typedef TMultiPhaseFilter MultiLevelSetType;
  typedef typename MultiLevelSetType::Pointer MultiLevelSetPointer;

  typedef ImageRegionIterator< ImageType > ImageIteratorType;
  typedef CastImageFilter< ImageType, OutputImageType > OutputCastType;
  typedef ImageFileWriter< InputImageType > WriterType;

  typedef ShapeLabelObject< ImagePixelType, ImageDimension >
    ShapeLabelObjectType;
  typedef LabelMap< ShapeLabelObjectType > ShapeLabelMapType;
  typedef typename ShapeLabelMapType::Pointer ShapeLabelMapPointer;
  typedef LabelImageToShapeLabelMapFilter< ImageType, ShapeLabelMapType >
    ShapeConverterType;
  typedef typename ShapeLabelMapType::LabelObjectContainerType
    LabelObjectContainerType;
  typedef typename LabelObjectContainerType::const_iterator LabelObjectIterator;

  typedef AtanHeavisideFunction< ImagePixelType, ImagePixelType > HeavisideFunctionType;

  typedef std::vector< unsigned int > VectorType;

  typedef Vector< float, ImageDimension > CentroidVectorType;
  typedef Statistics::ListSample< CentroidVectorType > SampleType;
  typedef Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer TreePointer;
  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef typename TreeType::Pointer KdTreePointer;

  itkGetConstMacro( Iterations, unsigned int );
  itkSetMacro( Iterations, unsigned int );
  itkGetConstMacro( RMSError, double );
  itkSetMacro( RMSError, double );
  itkGetConstMacro( Epsilon, double );
  itkSetMacro( Epsilon, double );
  itkGetConstMacro( CurvatureWeight, double );
  itkSetMacro( CurvatureWeight, double );
  itkGetConstMacro( AreaWeight, double );
  itkSetMacro( AreaWeight, double );
  itkGetConstMacro( Lambda1, double );
  itkSetMacro( Lambda1, double );
  itkGetConstMacro( Lambda2, double );
  itkSetMacro( Lambda2, double );
  itkGetConstMacro( OverlapPenaltyWeight, double );
  itkSetMacro( OverlapPenaltyWeight, double );
  itkGetConstMacro( LaplacianSmoothingWeight, double );
  itkSetMacro( LaplacianSmoothingWeight, double );
  itkGetConstMacro( VolumeMatchingWeight, double );
  itkSetMacro( VolumeMatchingWeight, double );
  itkGetConstMacro( Volume, double );
  itkSetMacro( Volume, double );
  itkGetConstMacro( LargestCellRadius, double );
  itkSetMacro( LargestCellRadius, double );
  itkGetConstMacro( FunctionCount, unsigned int );

  void SetFeatureImage( FeatureImageType *f )
  {
    m_Feature = f;
  }

  FeatureImageType *GetFeatureImage( void )
  {
    return m_Feature;
  }

  void SetSampling( float *sampling )
  {
    m_Sampling = sampling;
  }

protected:

  MultiLevelSetTracking();
  ~MultiLevelSetTracking(){}
  ImagePointer ComputeLevelSetFunction(
    ImagePointer input, unsigned int index );
  ImagePointer Resample ( ImagePointer iInput, ImageSpacingType iSpacing,
    ImageSizeType iSize, ImagePointType iOrigin );
  void GenerateData();

  typename ROIFilterType::Pointer m_ROIfilter;
  MultiLevelSetPointer m_LevelSetFilter;
  typename OutputCastType::Pointer m_CastOutput;

  unsigned int m_Iterations;
  double m_RMSError;
  double m_Epsilon;
  double m_CurvatureWeight;
  double m_AreaWeight;
  double m_Lambda1;
  double m_Lambda2;
  double m_OverlapPenaltyWeight;
  double m_LaplacianSmoothingWeight;
  double m_VolumeMatchingWeight;
  double m_Volume;
  unsigned int m_FunctionCount;
  double m_LargestCellRadius;
  float *m_Sampling;

  FeatureImageType *m_Feature;

private:

  MultiLevelSetTracking(Self&);   // intentionally not implemented
  void operator=(const Self&);   // intentionally not implemented
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiLevelSetTracking.txx"
#endif

#endif
