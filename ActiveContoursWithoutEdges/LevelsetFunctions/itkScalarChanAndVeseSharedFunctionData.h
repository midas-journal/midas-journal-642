#ifndef __itkScalarChanAndVeseSharedFunctionData_h_
#define __itkScalarChanAndVeseSharedFunctionData_h_

#include "itkLightObject.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkEuclideanDistance.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

template < class TInputImage, class TFeatureImage >
class ScalarChanAndVeseSharedFunctionData : public LightObject
{
public:

  typedef ScalarChanAndVeseSharedFunctionData Self;
  typedef LightObject Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int,
    TFeatureImage::ImageDimension );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  itkTypeMacro(ScalarChanAndVeseSharedFunctionData, LightObject);

  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::RegionType   InputRegionType;
  typedef typename InputImageType::SizeType     InputSizeType;
  typedef typename InputSizeType::SizeValueType InputSizeValueType;
  typedef typename InputImageType::SpacingType  InputSpacingType;
  typedef typename InputImageType::IndexType    InputIndexType;
  typedef typename InputIndexType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::PointType    InputPointType;

  typedef TFeatureImage FeatureImageType;
  typedef typename FeatureImageType::Pointer      FeatureImagePointer;
  typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;
  typedef typename FeatureImageType::PixelType    FeaturePixelType;
  typedef typename FeatureImageType::RegionType   FeatureRegionType;
  typedef typename FeatureImageType::SizeType     FeatureSizeType;
  typedef typename FeatureSizeType::SizeValueType FeatureSizeValueType;
  typedef typename FeatureImageType::SpacingType  FeatureSpacingType;
  typedef typename FeatureImageType::IndexType    FeatureIndexType;
  typedef typename FeatureImageType::PointType    FeaturePointType;

  typedef std::list< unsigned int > ListPixelType;
  typedef Image< ListPixelType, ImageDimension > ListImageType;
  typedef typename ListImageType::Pointer      ListImagePointer;
  typedef typename ListImageType::ConstPointer ListImageConstPointer;
  typedef typename ListImageType::RegionType   ListRegionType;
  typedef typename ListImageType::SizeType     ListSizeType;
  typedef typename ListSizeType::SizeValueType ListSizeValueType;
  typedef typename ListImageType::SpacingType  ListSpacingType;
  typedef typename ListImageType::IndexType    ListIndexType;
  typedef typename ListIndexType::IndexValueType ListIndexValueType;
  typedef typename ListImageType::PointType    ListPointType;
  typedef ImageRegionIteratorWithIndex< ListImageType > ListIteratorType;

  typedef Vector< float, ImageDimension > CentroidVectorType;
  typedef itk::Statistics::ListSample< CentroidVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer TreePointer;
  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef typename TreeType::Pointer KdTreePointer;

  void SetFunctionCount( const unsigned int& n )
  {
    FunctionCount = n;

    cVals.resize( n, 0. );
    cDens.resize( n, 0. );
    cNums.resize( n, 0. );

		cB.resize( n, 0. );
    cBDen.resize( n, 0. );
    cBNum.resize( n, 0. );

    //cBgrndNum = cBgrndDen = cBgrnd = 0;

    hVals.resize( n, 0 );
    start.resize( n );
    end.resize( n );
  }

  void CreateHVals( const unsigned int& j,
    const InputSpacingType& spacing,
    const InputPointType& origin,
    const InputRegionType& region )
  {
    hVals[j] = InputImageType::New();
    hVals[j]->SetRegions( region );
    hVals[j]->Allocate();
    hVals[j]->SetOrigin( origin );
    hVals[j]->SetSpacing( spacing );
    hVals[j]->FillBuffer( 0 );

    for( unsigned int i = 0; i < ImageDimension; i++ )
    {
      start[j][i] = static_cast< InputIndexValueType >( origin[i]/spacing[i] );
      end[j][i] = start[j][i] + static_cast< InputIndexValueType >(
        region.GetSize()[i] ) - 1;
    }
  }

  void SetKdTree( KdTreePointer kdtree )
  {
    this->m_KdTree = kdtree;
  }

  template< class TIndex >
  bool VerifyInsideRegion( const unsigned int& i, const TIndex& featureIndex )
  {
    typedef typename TIndex::IndexValueType TIndexValueType;
    for( unsigned int j = 0; j < ImageDimension; j++ )
    {
      if ((featureIndex[j] < static_cast< TIndexValueType >(start[i][j]) )
        || (featureIndex[j] > static_cast< TIndexValueType >(end[i][j])) )
        return false;
    }
    return true;
  }

  InputIndexType GetIndex( const unsigned int& j,
    const FeatureIndexType& featureIndex )
  {
    InputIndexType index;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      index[i] = featureIndex[i] -
        static_cast< InputIndexValueType >( start[j][i] );

		return index;
  }

  FeatureIndexType GetFeatureIndex( const unsigned int& j,
    const InputIndexType& inputIndex )
  {
    FeatureIndexType index;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      index[i] = inputIndex[i] +
        static_cast< InputIndexValueType >( start[j][i] );

    return index;
  }


  void AllocateListImage( FeatureRegionType region, FeatureSpacingType spacing )
  {
    lImage = ListImageType::New();
    lImage->SetRegions( region );
    lImage->Allocate();
    lImage->SetSpacing( spacing );
  }

  void PopulateListImage()
  {
    ListSpacingType spacing = lImage->GetSpacing();
    ListIteratorType lIt( lImage, lImage->GetLargestPossibleRegion() );

    if ( m_KdTree )
    {
      for(lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
      {
        ListIndexType ind = lIt.GetIndex();

        unsigned int numberOfNeighbors = 6;
        float queryPoint[ImageDimension];
        for( unsigned int i = 0; i < ImageDimension; i++ )
          queryPoint[i] = ind[i]*spacing[i];

        typename TreeType::InstanceIdentifierVectorType neighbors;
        this->m_KdTree->Search( queryPoint, numberOfNeighbors, neighbors ) ;

        ListPixelType L;
        for( unsigned int i = 0; i < numberOfNeighbors; i++ )
        {
          if ( VerifyInsideRegion( neighbors[i], ind ) )
            L.push_back( neighbors[i] );
        }
        lIt.Set( L );
      }
    }
    else
    {
      for(lIt.GoToBegin(); !lIt.IsAtEnd(); ++lIt )
      {
        ListIndexType ind = lIt.GetIndex();
        ListPixelType L;
        for( unsigned int i = 0; i < FunctionCount; i++ )
        {
          if ( VerifyInsideRegion( i, ind ) )
            L.push_back( i );
        }
        lIt.Set( L );
      }
    }
  }

  std::vector< double > cVals;
  std::vector< double > cNums;
  std::vector< double > cDens;
	std::vector< double > cB;
  std::vector< double > cBNum;
  std::vector< double > cBDen;

  unsigned int FunctionCount;
  std::vector< InputImagePointer >    hVals;
  std::vector< InputIndexType >  start;
  std::vector< InputIndexType >  end;
  ListImagePointer lImage;
  KdTreePointer m_KdTree;

protected:

  ScalarChanAndVeseSharedFunctionData()
  {
    m_KdTree = 0;
  }
  ~ScalarChanAndVeseSharedFunctionData(){}
};

} //end namespace itk

#endif
