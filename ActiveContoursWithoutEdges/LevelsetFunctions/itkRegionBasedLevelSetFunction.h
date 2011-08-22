#ifndef __itkRegionBasedLevelSetFunction_h_
#define __itkRegionBasedLevelSetFunction_h_

#include "itkFiniteDifferenceFunction.h"
#include "vnl/vnl_matrix_fixed.h"

#include "itkScalarChanAndVeseSharedFunctionData.h"
#include "itkSparseMultiphaseLevelSetImageFilter.h"
#include "itkDenseMultiphaseLevelSetImageFilter.h"
#include "itkRegularizedHeavisideFunctions.h"

namespace itk {

template < class TInput, // LevelSetImageType
  class TFeature, // FeatureImageType
  class TSharedData >
class ITK_EXPORT RegionBasedLevelSetFunction: public
FiniteDifferenceFunction< TInput >
{
public:
  /** Standard class typedefs. */
  typedef RegionBasedLevelSetFunction Self;
  typedef FiniteDifferenceFunction< TInput >
    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Run-time type information (and related methods) */
  itkTypeMacro( RegionBasedLevelSetFunction, FiniteDifferenceFunction );

  /** Extract some parameters from the superclass. */
  typedef double                                      TimeStepType;
  typedef typename Superclass::ImageType              ImageType;
  typedef typename Superclass::PixelType              PixelType;
  typedef PixelType ScalarValueType;
  typedef typename Superclass::RadiusType             RadiusType;
  typedef typename Superclass::NeighborhoodType       NeighborhoodType;
  typedef typename Superclass::NeighborhoodScalesType NeighborhoodScalesType;
  typedef typename Superclass::FloatOffsetType        FloatOffsetType;
  typedef FixedArray< ScalarValueType, ImageDimension > VectorType;

  struct GlobalDataStruct
  {
    ScalarValueType m_MaxAdvectionChange;
    ScalarValueType m_MaxPropagationChange;
    ScalarValueType m_MaxCurvatureChange;

    vnl_matrix_fixed<ScalarValueType,
      itkGetStaticConstMacro(ImageDimension),
      itkGetStaticConstMacro(ImageDimension)> m_dxy;

    ScalarValueType m_dx[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_dx_forward[itkGetStaticConstMacro(ImageDimension)];
    ScalarValueType m_dx_backward[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_GradMagSqr;
  };

  /* This structure is derived from LevelSetFunction and stores intermediate
  values for computing time step sizes */
  struct ACGlobalDataStruct : public GlobalDataStruct
  {
    ScalarValueType m_MaxAdvectionChange;
    ScalarValueType m_MaxPropagationChange;
    ScalarValueType m_MaxCurvatureChange;

    vnl_matrix_fixed<ScalarValueType,
      itkGetStaticConstMacro(ImageDimension),
      itkGetStaticConstMacro(ImageDimension)> m_dxy;

    ScalarValueType m_dx[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_dx_forward[itkGetStaticConstMacro(ImageDimension)];
    ScalarValueType m_dx_backward[itkGetStaticConstMacro(ImageDimension)];

    ScalarValueType m_GradMagSqr;
    ScalarValueType m_MaxGlobalChange;
  };


  typedef TInput InputImageType;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename InputImageType::IndexType InputIndexType;
  typedef typename InputImageType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename InputImageType::SizeValueType InputSizeValueType;
  typedef typename InputImageType::RegionType InputRegionType;
  typedef typename InputImageType::PointType InputPointType;

  typedef TFeature FeatureImageType;
  typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;
  typedef typename FeatureImageType::PixelType FeaturePixelType;
  typedef typename FeatureImageType::IndexType FeatureIndexType;
  typedef typename FeatureImageType::OffsetType FeatureOffsetType;

  typedef TSharedData SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

  typedef ImageRegionIteratorWithIndex< InputImageType >
    ImageIteratorType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType >
    ConstImageIteratorType;
  typedef ImageRegionIteratorWithIndex< FeatureImageType >
    FeatureImageIteratorType;
  typedef ImageRegionConstIterator< FeatureImageType >
    ConstFeatureIteratorType;

  typedef HeavisideFunctionBase< InputPixelType, InputPixelType > HeavisideFunctionType;

  void SetDomainFunction( HeavisideFunctionType* f )
  {
    m_DomainFunction = f;
  }

  virtual void Initialize(const RadiusType &r)
  {
    this->SetRadius(r);

    // Dummy neighborhood.
    NeighborhoodType it;
    it.SetRadius( r );

    // Find the center index of the neighborhood.
    m_Center =  it.Size() / 2;

    // Get the stride length for each axis.
    for(unsigned int i = 0; i < ImageDimension; i++)
      {  m_xStride[i] = it.GetStride(i); }
  }



  void SetSharedData( SharedDataPointer sharedDataIn )
  {
    sharedData = sharedDataIn;
  }

  void UpdateSharedData( bool forceUpdate );

  void *GetGlobalDataPointer() const
  {
    ACGlobalDataStruct *ans = new ACGlobalDataStruct();

    ScalarValueType null_value = NumericTraits<ScalarValueType>::Zero;

    ans->m_MaxCurvatureChange   = null_value;
	  ans->m_MaxGlobalChange      = null_value;
    return ans;
  }

  TimeStepType ComputeGlobalTimeStep(void *GlobalData) const;

  /** Compute the equation value. */
  PixelType ComputeUpdate(const NeighborhoodType &neighborhood,
    void *globalData, const FloatOffsetType& = FloatOffsetType(0.0));

  void SetInitialImage(InputImageType *f)
  {
    m_InitialImage = f;
  }

  virtual const FeatureImageType *GetFeatureImage() const
    { return m_FeatureImage.GetPointer(); }
  virtual void SetFeatureImage(const FeatureImageType *f)
    {    m_FeatureImage = f;  }

  /** Nu. Area regularization values */
  void SetAreaWeight( const ScalarValueType& nu)
    { this->m_AreaWeight = nu; }
  ScalarValueType GetAreaWeight() const
    { return this->m_AreaWeight; }

  /** Lambda1. Internal intensity difference weight */
  void SetLambda1( const ScalarValueType& lambda1 )
    { this->m_Lambda1 = lambda1; }
  ScalarValueType GetLambda1() const
    { return this->m_Lambda1; }

  /** Lambda2. External intensity difference weight */
  void SetLambda2( const ScalarValueType& lambda2 )
    { this->m_Lambda2 = lambda2; }
  ScalarValueType GetLambda2() const
    { return this->m_Lambda2; }

  /** Gamma. Overlap penalty */
  void SetOverlapPenaltyWeight( const ScalarValueType& gamma )
    { this->m_OverlapPenaltyWeight = gamma; }
  ScalarValueType GetOverlapPenaltyWeight() const
    { return this->m_OverlapPenaltyWeight; }

  /** Gamma. Scales all curvature weight values */
  virtual void SetCurvatureWeight(const ScalarValueType c)
    { m_CurvatureWeight = c; }
  ScalarValueType GetCurvatureWeight() const
    { return m_CurvatureWeight; }

  /** Weight of the laplacian smoothing term */
  void SetLaplacianSmoothingWeight(const ScalarValueType c)
    { m_LaplacianSmoothingWeight = c; }
  ScalarValueType GetLaplacianSmoothingWeight() const
    { return m_LaplacianSmoothingWeight; }

  /** Volume matching weight.  */
  void SetVolumeMatchingWeight( const ScalarValueType& tau )
    { this->m_VolumeMatchingWeight = tau; }
  ScalarValueType GetVolumeMatchingWeight() const
    { return this->m_VolumeMatchingWeight; }

  /** Volume.  */
  void SetVolume( const ScalarValueType& volume )
    { this->m_Volume = volume; }
  ScalarValueType GetVolume() const
    { return this->m_Volume; }

  /** Set function id.  */
  void SetFunctionId( const unsigned int& iFid )
    { this->m_FunctionId = iFid; }

  virtual void ReleaseGlobalDataPointer(void *GlobalData) const
  { delete (GlobalDataStruct *) GlobalData; }

  virtual ScalarValueType ComputeCurvatureTerm(const NeighborhoodType &,
    const FloatOffsetType &, GlobalDataStruct *gd = 0 );

  /** Laplacian smoothing speed.  Can be used to spatially modify the
    effects of laplacian smoothing of the level set function */
  virtual ScalarValueType LaplacianSmoothingSpeed(
    const NeighborhoodType &,
    const FloatOffsetType &, GlobalDataStruct * = 0) const
    { return NumericTraits<ScalarValueType>::One; }

  /** Curvature speed.  Can be used to spatially modify the effects of
    curvature . The default implementation returns one. */
  virtual ScalarValueType CurvatureSpeed(const NeighborhoodType &,
                                         const FloatOffsetType &, GlobalDataStruct * = 0
                                         ) const
    { return NumericTraits<ScalarValueType>::One; }

protected:

  RegionBasedLevelSetFunction();
  virtual ~RegionBasedLevelSetFunction() {}

  /** The initial level set image */
  InputImageConstPointer m_InitialImage;

  /** The feature image */
  FeatureImageConstPointer m_FeatureImage;

  bool updatedC;
  bool updatedH;

  SharedDataPointer sharedData;
  HeavisideFunctionType* m_DomainFunction;

  /* Area regularizer term in CV formulation, what about lambda1 and lambda2?*/
  ScalarValueType m_AreaWeight;
  ScalarValueType m_Lambda1;
  ScalarValueType m_Lambda2;
  ScalarValueType m_OverlapPenaltyWeight;
	ScalarValueType m_VolumeMatchingWeight;
	ScalarValueType m_Volume;
  unsigned int m_FunctionId;

  std::slice x_slice[itkGetStaticConstMacro(ImageDimension)];
  ::size_t m_Center;
  ::size_t m_xStride[itkGetStaticConstMacro(ImageDimension)];

  static double m_WaveDT;
  static double m_DT;

  ScalarValueType m_CurvatureWeight;
  ScalarValueType m_LaplacianSmoothingWeight;


  void ComputeHImage();
  ScalarValueType computeRegularizationTerms( /* fill with adequate parameters */);

  ScalarValueType computeGlobalTerm(
    const ScalarValueType& imagePixel,
    const InputIndexType& inputIndex );

  virtual ScalarValueType computeInternalTerm(const FeaturePixelType& iValue,
    const FeatureIndexType& iIdx, const unsigned int& fId ) = 0;

  virtual ScalarValueType computeExternalTerm(const FeaturePixelType& iValue,
    const FeatureIndexType& iIdx, const unsigned int& pr ) = 0;

  virtual void computeOverlapParameters( const FeatureIndexType featIndex,
    unsigned int& s, unsigned int& pr ) = 0;

  virtual ScalarValueType computeOverlapTerm( const unsigned int& s )
  { return this->m_OverlapPenaltyWeight * s; }

  virtual void ComputeParameters() = 0;

  virtual void SpecialProcessing(){}

private:
  RegionBasedLevelSetFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRegionBasedLevelSetFunction.txx"
#endif

#endif
