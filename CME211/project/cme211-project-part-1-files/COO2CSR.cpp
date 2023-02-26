#include <algorithm>
#include <vector>

#include "COO2CSR.hpp"

/* Bubble sorting function for sparse matrix format conversion,
   used to sort the entries in one row of the matrix. */

void SortRow(std::vector<int> &row_idx, std::vector<int> &col_idx,
             std::vector<double> &a, unsigned int start, unsigned int end, unsigned int nz) {
  for (unsigned int i = end - 1; i > start; i--) {
    for (unsigned int j = start; j < i; j++) {
      if (col_idx[j] > col_idx[j + 1]) {
        /* Swap the value and the column index */
        double dt = a[j];
        a[j] = a[j + 1];
        a[j + 1] = dt;

        int it = col_idx[j];
        col_idx[j] = col_idx[j + 1];
        col_idx[j + 1] = it;
      }
    }
  }

  /* Accumulate duplicate values and adjust vectors */
  for (unsigned int j = start; j < end - 1; j++) {
    if (col_idx[j] == col_idx[j + 1]) {
      a[j] += a[j + 1];
      for (unsigned int i = j + 1; i < nz - 1; i++) {
        a[i] = a[i + 1];
        col_idx[i] = col_idx[i + 1];
        if (row_idx[i] > 0)
          row_idx[i]--;
      }
      col_idx[nz - 1] = -1;
      a[nz - 1] = 0;
      end--;
      j--;
    }
  }
}

/* In place conversion of square matrix from COO to CSR format */

void COO2CSR(std::vector<double> &val, std::vector<int> &i_idx,
             std::vector<int> &j_idx) {
  unsigned int nrows = (unsigned int) *std::max_element(i_idx.begin(), i_idx.end()) + 1;
  unsigned int nnonzeros = (unsigned int)val.size();

  unsigned int *row_start = new unsigned int[nrows + 1];
  for (unsigned int i = 0; i <= nrows; i++) {
    row_start[i] = 0;
  }

  /* Determine row lengths */

  for (unsigned int i = 0; i < nnonzeros; i++)
    row_start[(unsigned int)i_idx[i] + 1]++;
  for (unsigned int i = 0; i < nrows; i++)
    row_start[i + 1] += row_start[i];

  for (unsigned int init = 0; init < nnonzeros;) {
    double dt = val[init];
    int i = i_idx[init];
    int j = j_idx[init];
    i_idx[init] = -1;

    while (true) {
      unsigned int i_pos = row_start[i];
      double val_next = val[i_pos];
      int i_next = i_idx[i_pos];
      int j_next = j_idx[i_pos];

      val[i_pos] = dt;
      j_idx[i_pos] = (int) j;
      i_idx[i_pos] = -1;
      row_start[i]++;

      if (i_next < 0)
        break;

      dt = val_next;
      i = i_next;
      j = j_next;
    }
    init++;

    while ((init < nnonzeros) and (i_idx[init] < 0)) {
      init++;
    }
  }

  /* Copy row pointer */

  for (unsigned int i = 0; i < nrows; i++) {
    i_idx[i + 1] = (int) row_start[i];
  }
  i_idx[0] = 0;

  /* Sort each row */

  for (unsigned int i = 0; i < nrows; i++) {
    SortRow(i_idx, j_idx, val, (unsigned int) i_idx[i], (unsigned int) i_idx[i + 1], nnonzeros);
  }

  i_idx.resize(nrows + 1);

  delete[] row_start;
}
