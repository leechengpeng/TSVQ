#include <iostream>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../../TSVQ.hpp"

int main()
{
	cv::Mat myImage = cv::imread("waterfall.jpg");

	auto ImageWidth = myImage.cols;
	auto ImageHeight = myImage.rows;
	auto PixelsSize = ImageHeight * ImageWidth;

	unsigned int Dimension = 0;
	for (int i=3; i<PixelsSize; ++i)
	{
		if (PixelsSize % i)
		{
			Dimension = i;
			break;
		}
	}

	std::vector<uchar*> ImageSet;
	for (int rows=0; rows<ImageHeight; ++rows)
	{
		for (int cols=0; cols<ImageWidth; ++cols)
		{
			ImageSet.push_back(myImage.ptr(rows) + cols * 3);
		}
	}

	CTSVQ<uchar> TSVQ;
	TSVQ.quantizeVectors(ImageSet, 3, 2);

	for (int rows=0; rows<ImageHeight; ++rows)
	{
		for (int cols=0; cols<ImageWidth; ++cols)
		{
			auto pRows = myImage.ptr(rows);
			auto pPixel = TSVQ.findCodeVectors(pRows + cols * 3);
			for (int i=0; i<3; ++i)
			{
				pRows[cols*3 + i] = *(pPixel + i);
			}
		}
	}

	cv::imshow("test", myImage);
	cv::waitKey();

	return 0;
}