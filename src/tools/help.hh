/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef __HELP_HH
#define __HELP_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

/**
 * class Float1D is a proxy class that can represent, for example, 
 * a column or row vector of a Float2D array, where row (sub-)arrays 
 * are stored with a respective stride. 
 * Besides constructor/deconstructor, the class provides overloading of 
 * the []-operator, such that elements can be accessed as v[i] 
 * (independent of the stride).
 * The class will never allocate separate memory for the vectors, 
 * but point to the interior data structure of Float2D (or other "host" 
 * data structures).
 */ 
class Float1D
{
  public:
	Float1D(float* _elem, int _rows, int _stride = 1) 
	: rows(_rows),stride(_stride),elem(_elem)
	{
	}

	Float1D(int _rows, int _stride = 1)
	: rows(_rows),stride(_stride)
	{
		elem = new float[rows*stride]();
	}

	~Float1D()
	{
	}

	inline float& operator[](int i) { 
		return elem[i*stride]; 
	}

	inline const float& operator[](int i) const {
		return elem[i*stride]; 
	}

	inline float* elemVector() {
		return elem;
	}

        inline int getSize() const { return rows; }; 

  private:
    int rows;
    int stride;
    float* elem;
};

/**
 * class Float2D is a very basic helper class to deal with 2D float arrays:
 * indices represent columns (1st index, "horizontal"/x-coordinate) and 
 * rows (2nd index, "vertical"/y-coordinate) of a 2D grid;
 * values are sequentially ordered in memory using "column major" order.
 * Besides constructor/deconstructor, the class provides overloading of 
 * the []-operator, such that elements can be accessed as a[i][j]. 
 */ 
class Float2D
{
public:
	/**
	 * Constructor
	 * @param _cols	number of columns (i.e., elements in horizontal direction)
	 * @param _rows rumber of rows (i.e., elements in vertical directions)
	 */
	Float2D(int _cols, int _rows, int _cpadd=0, int  _rpadd=0)
    : rows(_rows), cols(_cols), cpadd(_cpadd), rpadd(_rpadd), isProxy(false)
	{
		elem = new float[(rows+rpadd)*(cols+cpadd)];
		for(int i=0; i < cols + cpadd; i++)
			for(int j=0; j < rows + rpadd; j++)
				elem[i*(rows+rpadd)+j] = 0;
	}

	Float2D(float* _elem, int _cols, int _rows, int _cpadd=0, int  _rpadd=0)
	: rows(_rows), cols(_cols), elem(_elem), cpadd(_cpadd), rpadd(_rpadd), isProxy(true) {}

    Float2D(const Float2D& obj)
    : rows(obj.rows), cols(obj.cols), cpadd(obj.cpadd), rpadd(obj.rpadd), isProxy(obj.isProxy)
    {
		if (isProxy)
			elem = obj.elem;
		else {
			elem = new float[(rows+rpadd)*(cols+cpadd)];
			for(int i=0; i < cols + cpadd; i++)
				for(int j=0; j < rows + rpadd; j++)
					elem[i*(rows+rpadd)+j] = obj.elem[i*(rows+rpadd)+j];
		}
    }

	~Float2D()
	{
		if (!isProxy)
			delete[] elem;
	}

	inline float* operator[](int i) {
		return (elem + ((rows+rpadd) * i));
	}

	inline float const* operator[](int i) const {
		return (elem + ((rows+rpadd) * i));
	}

	inline float* elemVector() {
		return elem;
	}

    inline int getRows() const { return rows; };
    inline int getCols() const { return cols; };

	inline Float1D getColProxy(int i) {
		// subarray elem[i][*]:
        // starting at elem[i][0] with rows elements and unit stride
		return Float1D(elem + ((rows+rpadd) * i), rows);
	};

	inline Float1D getRowProxy(int j) {
		// subarray elem[*][j]
        // starting at elem[0][j] with cols elements and stride rows
		return Float1D(elem + j, cols, rows+rpadd);
	};

	inline Float2D* getBlockProxy(int c, int r, int nc, int nr) {
		// submatrix starting at elem[c][r] with nr rows and nc columns
//		Float2D proxy(nc, nr);
//		for (int i=0; i<nc; i++)
//			for (int j=0; j<nr; j++)
//				proxy[i][j] = elem[(i+c)*rows+j+r];
//		return proxy;
		return new Float2D(elem + c*(rows+rpadd) + r, nc, nr, cols+cpadd-nc, rows+rpadd-nr);
	}

  private:
    int rows;
    int cols;
    float* elem;
    // padding for rows and columns
	int cpadd;
	int rpadd;
	bool isProxy;
};

/*
 * returns the smallest if all are positive,
           the largest  if all are negative,
           zero         if they are not all positive or all negative
 */
inline float minmod(float a, float b) {
	if (a>0 && b>0)
		return std::min(a,b);
	if (a<0 && b<0)
		return std::max(a,b);
	return 0.;
}

//-------- Methods for Visualistion of Results --------

/**
 * generate output filenames for the single-SWE_Block version
 * (for serial and OpenMP-parallelised versions that use only a 
 *  single SWE_Block - one output file is generated per checkpoint)
 */
inline std::string generateFileName(std::string baseName, int timeStep) {

	std::ostringstream FileName;
	FileName << baseName <<timeStep<<".vtk";
	return FileName.str();
};

/**
 * Generates an output file name for a multiple SWE_Block version based on the ordering of the blocks.
 *
 * @param i_baseName base name of the output.
 * @param i_blockPositionX position of the SWE_Block in x-direction.
 * @param i_blockPositionY position of the SWE_Block in y-direction.
 * @param i_fileExtension file extension of the output file.
 * @return
 */
inline std::string generateFileName( std::string i_baseName,
                                     int i_blockPositionX, int i_blockPositionY,
                                     std::string i_fileExtension=".nc" ) {

  std::ostringstream l_fileName;

  l_fileName << i_baseName << "_" << i_blockPositionX << i_blockPositionY << i_fileExtension;
  return l_fileName.str();
};

/**
 * generate output filename for the multiple-SWE_Block version
 * (for serial and parallel (OpenMP and MPI) versions that use 
 *  multiple SWE_Blocks - for each block, one output file is 
 *  generated per checkpoint)
 */
inline std::string generateFileName(std::string baseName, int timeStep, int block_X, int block_Y, std::string i_fileExtension=".vts") {

	std::ostringstream FileName;
	FileName << baseName <<"_"<< block_X<<"_"<<block_Y<<"_"<<timeStep<<i_fileExtension;
	return FileName.str();
};

/**
 * generate output filename for the ParaView-Container-File
 * (to visualize multiple SWE_Blocks per checkpoint)
 */
inline std::string generateContainerFileName(std::string baseName, int timeStep) {

	std::ostringstream FileName;
	FileName << baseName<<"_"<<timeStep<<".pvts";
	return FileName.str();
};


#endif

