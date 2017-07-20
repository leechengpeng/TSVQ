#include <iostream>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../../TSVQ.hpp"

int main()
{
	cv::Mat Image = cv::imread("waterfall.jpg");

	size_t ImageWidth = Image.cols, ImageHeight = Image.rows;
	size_t PixelsSize = ImageHeight * ImageWidth;

	std::vector<uchar*> ImageSet;
	for (int rows=0; rows<ImageHeight; ++rows)
	{
		for (int cols=0; cols<ImageWidth; ++cols)
		{
			ImageSet.push_back(Image.ptr(rows) + cols * 3);
		}
	}

	CTSVQ<uchar> TSVQ;
	TSVQ.quantizeVectors(ImageSet, Image.channels(), 16);

	for (int rows=0; rows<ImageHeight; ++rows)
	{
		for (int cols=0; cols<ImageWidth; ++cols)
		{
			auto pRows = Image.ptr(rows);
			auto pPixel = TSVQ.findCodeVectors(pRows + cols * 3);
			for (int i=0; i<3; ++i)
			{
				pRows[cols*3 + i] = *(pPixel + i);
			}
		}
	}
	cv::imwrite("VQ.jpg", Image);

	return 0;
}