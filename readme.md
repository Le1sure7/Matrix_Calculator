IdentityMatrix(const int Size) Make a identity matrix whose size is Size;
/name num Name
/inverse num/Name
/transpose num/Name
/determinant num/Name
/rank num/Name
/read
/cal num/Name + num/Name
/cal num/Name - num/Name
/cal num/Name @ num/Name
/cal num/Name ^ c
/cal c * num/Name
/cal num/Name * c
/cal rowOP     
    // Rows are numbered from 1.
    // type 1: swap row1 and row2.
    // type 2: multiply row1 by factor (nonzero).
    // type 3: row1 += factor * row2 (distinct rows).
    RowTransform(const int type,const int row1,const int row2 = 0,const double factor = 1.0)