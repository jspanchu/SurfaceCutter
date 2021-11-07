/**
MIT License

Copyright (c) 2021 Jaswant Sai Panchumarti

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#define WRITE_TSC_TEST_DATA 0

#include "tscTriSurfaceCutter.h"
#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkDataArray.h>
#include <vtkDataArrayAccessor.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataReader.h>
#if WRITE_TSC_TEST_DATA
#include <vtkXMLMultiBlockDataWriter.h>
#endif

#include <array>
#include <sstream>
#include <string>

using dispatchRR =
  vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;
using dispatchAA =
  vtkArrayDispatch::Dispatch2ByValueType<vtkArrayDispatch::AllTypes, vtkArrayDispatch::AllTypes>;

struct DataArrComparator
{
  template <typename ScalarsArr1T, typename = vtk::detail::EnableIfVtkDataArray<ScalarsArr1T>,
    typename ScalarsArr2T, typename = vtk::detail::EnableIfVtkDataArray<ScalarsArr2T>>
  void operator()(ScalarsArr1T* scalarsArr1, ScalarsArr2T* scalarsArr2, bool& success)
  {
    success = false;

    vtkIdType nvals(0);
    if (!(nvals = scalarsArr1->GetNumberOfValues()))
      return;

    if (nvals != scalarsArr2->GetNumberOfValues())
      return;

    auto scalars1 = vtk::DataArrayValueRange(scalarsArr1);
    auto scalars2 = vtk::DataArrayValueRange(scalarsArr2);

    using Scalars1T = vtk::GetAPIType<ScalarsArr1T>;
    using Scalars2T = vtk::GetAPIType<ScalarsArr2T>;

    for (vtkIdType val = 0; val < nvals; ++val)
    {
      if (sizeof(Scalars1T) <= sizeof(Scalars2T))
      {
        constexpr auto eps = std::numeric_limits<Scalars1T>::epsilon();
        if (scalars1[val] - scalars2[val] > eps)
          return;
      }
      else
      {
        constexpr auto eps = std::numeric_limits<Scalars2T>::epsilon();
        if (scalars1[val] - scalars2[val] > eps)
          return;
      }
    }

    success = true;
  }
};

int compare(vtkSmartPointer<vtkPolyData> pdata1, vtkSmartPointer<vtkPolyData> pdata2)
{

  if (!(pdata1->GetPoints() || pdata2->GetPoints()))
  {
    std::cerr << "Unexpected error! Empty point set." << std::endl;
    return EXIT_FAILURE;
  }

  vtkDataArray* points1 = pdata1->GetPoints()->GetData();
  vtkDataArray* points2 = pdata2->GetPoints()->GetData();

  bool success(false);
  auto comparator = DataArrComparator();
  if (!dispatchRR::Execute(points1, points2, comparator, success))
    comparator(points1, points2, success);

  if (!success)
    return EXIT_FAILURE;

  int numArrs1 = pdata1->GetPointData()->GetNumberOfArrays();
  int numArrs2 = pdata2->GetPointData()->GetNumberOfArrays();
  if (numArrs1 != numArrs2)
    return EXIT_FAILURE;
  for (int iArr = 0; iArr < numArrs1; ++iArr)
  {
    success = false;
    vtkDataArray* arr1 = pdata1->GetPointData()->GetArray(iArr);
    vtkDataArray* arr2 = pdata2->GetPointData()->GetArray(iArr);

    if (!dispatchAA::Execute(arr1, arr2, comparator, success))
      comparator(arr1, arr2, success);

    if (!success)
      return EXIT_FAILURE;
  }

  numArrs1 = pdata1->GetCellData()->GetNumberOfArrays();
  numArrs2 = pdata2->GetCellData()->GetNumberOfArrays();
  if (numArrs1 != numArrs2)
    return EXIT_FAILURE;
  for (int iArr = 0; iArr < numArrs1; ++iArr)
  {
    success = false;
    vtkDataArray* arr1 = pdata1->GetCellData()->GetArray(iArr);
    vtkDataArray* arr2 = pdata2->GetCellData()->GetArray(iArr);

    if (!dispatchAA::Execute(arr1, arr2, comparator, success))
      comparator(arr1, arr2, success);

    if (!success)
      return EXIT_FAILURE;
  }

  const vtkIdType& numTris = pdata1->GetNumberOfPolys();
  if (numTris != pdata2->GetNumberOfPolys())
    return EXIT_FAILURE;

  const vtkIdType& numLines = pdata1->GetNumberOfLines();
  if (numLines != pdata2->GetNumberOfLines())
    return EXIT_FAILURE;

  auto iter = vtk::TakeSmartPointer(pdata1->NewCellIterator());
  auto testIter = vtk::TakeSmartPointer(pdata2->NewCellIterator());
  for (iter->InitTraversal(), testIter->InitTraversal();
       !iter->IsDoneWithTraversal() && !testIter->IsDoneWithTraversal();
       iter->GoToNextCell(), testIter->GoToNextCell())
  {
    // same cell type
    const int& cellType = iter->GetCellType();
    const int& testCellType = testIter->GetCellType();
    if (cellType != testCellType)
      return EXIT_FAILURE;

    // same pt ids.
    vtkSmartPointer<vtkIdList> ptIds = iter->GetPointIds();
    vtkSmartPointer<vtkIdList> testPtIds = testIter->GetPointIds();
    const vtkIdType& numPtIds = ptIds->GetNumberOfIds();
    const vtkIdType& numTestPtIds = testPtIds->GetNumberOfIds();
    if (numPtIds != numTestPtIds)
      return EXIT_FAILURE;
    for (vtkIdType iPtId = 0; iPtId < numPtIds; ++iPtId)
    {
      if (ptIds->GetId(iPtId) != testPtIds->GetId(iPtId))
      {
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}

int main()
{
  vtkNew<tscTriSurfaceCutter> surfCutter;
  vtkNew<vtkPolyData> triangle;
  vtkNew<vtkPoints> tpoints;
  vtkNew<vtkCellArray> tCell;
  vtkNew<vtkXMLMultiBlockDataReader> reader;
#if WRITE_TSC_TEST_DATA
  vtkNew<vtkXMLMultiBlockDataWriter> writer;
#endif
  vtkNew<vtkMultiBlockDataSet> inOutTrueBlks, inOutFalseBlks;

  triangle->SetPoints(tpoints);
  triangle->SetPolys(tCell);

  // init triangle
  tpoints->InsertNextPoint(-1.0, -1.0, 0.0);
  tpoints->InsertNextPoint(1.0, -1.0, 0.0);
  tpoints->InsertNextPoint(1.0, 1.0, 0.0);
  tCell->InsertNextCell({ 0, 1, 2 });

  surfCutter->SetInputData(0, triangle);

  // data for loops
  std::vector<std::vector<std::vector<double>>> testPoints = {
    { { -1.01951801776886, -0.2836509943008423, 0.0 }, { 1.0, -1.0, 0.0 },
      { 0.3121359944343567, 1.033759951591492, 0.0 } }, // l0
    { { -1.462092995643616, -1.623507976531982, 0.0 }, { 1.0, -1.0, 0.0 },
      { -0.5603539943695068, 0.7119830250740051, 0.0 } }, // l1
    { { -1.18078601360321, -0.6734970211982727, 0.0 },
      { 1.246237993240356, -1.113623976707458, 0.0 },
      { 0.5600630044937134, 1.047474980354309, 0.0 } }, // l2
    { { -1.463358998298645, -0.5595560073852539, 0.0 },
      { 1.770364999771118, -0.001563000027090311, 0.0 },
      { -1.454408049583435, 0.591713011264801, 0.0 } }, // l3
    { { 1.770364999771118, -0.001563000027090311, 0.0 },
      { -2.11070704460144, 0.541579008102417, 0.0 },
      { 0.8609480261802673, -1.524399995803833, 0.0 },
      { -0.4971419870853424, -0.04454400017857552, 0.0 } },      // l4
    { { 1.0, 0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 0.0, 0.0, 0.0 } } // l5
  };
  std::vector<std::vector<vtkIdType>> testPolys = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
    { 0, 1, 2 }, { 0, 1, 2, 3 }, { 0, 1, 2 } };

  reader->SetFileName("Data/InOutTrue.vtm");
  reader->Update();
  inOutTrueBlks->ShallowCopy(reader->GetOutput());
  for (unsigned short i = 0; i < 6; ++i)
  {
    vtkNew<vtkPolyData> loop;
    vtkNew<vtkPoints> lPoints;
    vtkNew<vtkCellArray> lCell;
    for (const auto& point : testPoints[i])
    {
      lPoints->InsertNextPoint(point.data());
    }
    lCell->InsertNextCell(testPolys[i].size(), testPolys[i].data());
    loop->SetPoints(lPoints);
    loop->SetPolys(lCell);
    surfCutter->SetLoopsData(loop);
    surfCutter->Update();

    vtkSmartPointer<vtkPolyData> cut = surfCutter->GetOutput();
#if WRITE_TSC_TEST_DATA
    vtkNew<vtkPolyData> test;
    test->DeepCopy(cut);
    inOutTrueBlks->SetBlock(i, test);
#else
    vtkSmartPointer<vtkPolyData> test = vtkPolyData::SafeDownCast(inOutTrueBlks->GetBlock(i));
#endif

    if (compare(cut, test) == EXIT_FAILURE)
    {
      std::cerr << "InsideOut: True | Test " << i << ": Failed\n";
      return EXIT_FAILURE;
    }

    std::cout << "InsideOut: True | Test " << i << ": Passed\n";
  }
#if WRITE_TSC_TEST_DATA
  writer->SetInputData(inOutTrueBlks);
  writer->SetFileName("InOutTrue.vtm");
  writer->Write();
#endif
  surfCutter->SetInsideOut(false);

  reader->SetFileName("Data/InOutFalse.vtm");
  reader->Update();
  inOutFalseBlks->ShallowCopy(reader->GetOutput());
  for (unsigned short i = 0; i < 6; ++i)
  {
    vtkNew<vtkPolyData> loop;
    vtkNew<vtkPoints> lPoints;
    vtkNew<vtkCellArray> lCell;
    for (const auto& point : testPoints[i])
    {
      lPoints->InsertNextPoint(point.data());
    }
    lCell->InsertNextCell(testPolys[i].size(), testPolys[i].data());
    loop->SetPoints(lPoints);
    loop->SetPolys(lCell);
    surfCutter->SetLoopsData(loop);
    surfCutter->Update();

    vtkSmartPointer<vtkPolyData> cut = surfCutter->GetOutput();
#if WRITE_TSC_TEST_DATA
    vtkNew<vtkPolyData> test;
    test->DeepCopy(cut);
    inOutFalseBlks->SetBlock(i, test);
#else
    vtkSmartPointer<vtkPolyData> test = vtkPolyData::SafeDownCast(inOutFalseBlks->GetBlock(i));
#endif

    if (compare(cut, test) == EXIT_FAILURE)
    {
      std::cerr << "InsideOut: False | Test " << i << ": Failed\n";
      return EXIT_FAILURE;
    }

    std::cout << "InsideOut: False| Test " << i << ": Passed\n";
  }
#if WRITE_TSC_TEST_DATA
  writer->SetInputData(inOutFalseBlks);
  writer->SetFileName("InOutFalse.vtm");
  writer->Write();
#endif

  return EXIT_SUCCESS;
}