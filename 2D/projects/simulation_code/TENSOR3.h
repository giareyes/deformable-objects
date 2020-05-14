#include "EXTRAFUNCTIONS.h"
#include <vector>

using namespace std;

class TENSOR3
{
public:
  TENSOR3();
  TENSOR3(int rows, int cols, int slabs);
  TENSOR3(const std::vector<MATRIX>& slabs);

  int rows() const { return _rows; };
  int cols() const { return _cols; };
  int slabs() const { return _slabs; };


  MATRIX modeThreeProduct(const VECTOR& x);
  TENSOR3 modeThreeProduct(const MATRIX& x);

  MATRIX modeTwoProduct(const VECTOR& x);
  TENSOR3 modeTwoProduct(const MATRIX& x);

  TENSOR3 modeOneProduct(const MATRIX& x);
  TENSOR3& operator+=(const TENSOR3& m);
  TENSOR3& operator*=(const Real& scalar);

  void toString();

  void clear();

  vector<MATRIX> _tensor;
protected:
  int _rows;
  int _cols;
  int _slabs;
};
