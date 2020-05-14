#include "TENSOR3.h"
// #include <vector>

class TENSOR4
{
public:
  TENSOR4();
  TENSOR4(int rows, int cols, int slab_rows, int slab_cols);
  TENSOR4(const vector<TENSOR3>& slabs);

  int rows() const { return _rows; };
  int cols() const { return _cols; };
  int slab_rows() const { return _slab_rows; };
  int slab_cols() const { return _slab_cols; };

  void toString();
  void clear();

  TENSOR3 modeFourProduct(const VECTOR& x);
  TENSOR4 modeFourProduct(const MATRIX& x);

  TENSOR4 modeThreeProduct(const MATRIX& x);
  TENSOR3 modeThreeProduct(const VECTOR& x);

  TENSOR4 modeTwoProduct(const MATRIX& x);
  TENSOR3 modeTwoProduct(const VECTOR& x);

  TENSOR4 modeOneProduct(const MATRIX& x);

  TENSOR4& operator+=(const TENSOR4& m);
  TENSOR4& operator*=(const Real& scalar);

  vector<TENSOR3> _tensor;
protected:
  int _rows;
  int _cols;
  int _slab_rows;
  int _slab_cols;

};
