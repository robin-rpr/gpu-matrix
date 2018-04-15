import LuDecomposition from './dc/lu';
import SvDecomposition from './dc/svd';
import rescale from 'ml-array-rescale';
import {
    checkRowVector, checkRowIndex, checkColumnIndex, checkColumnVector,
    checkRange, checkIndices,
    sumByRow, sumByColumn, sumAll
} from './util';
import MatrixTransposeView from './views/transpose';
import MatrixRowView from './views/row';
import MatrixSubView from './views/sub';
import MatrixSelectionView from './views/selection';
import MatrixRowSelectionView from './views/rowSelection';
import MatrixColumnSelectionView from './views/columnSelection';
import MatrixColumnView from './views/column';
import MatrixFlipRowView from './views/flipRow';
import MatrixFlipColumnView from './views/flipColumn';

export default function AbstractMatrix(superCtor) {
    if (superCtor === undefined) superCtor = Object;

    /**
     * Real matrix
     * @class Matrix
     * @param {number|Array|Matrix} nRows - Number of rows of the new matrix,
     * 2D array containing the data or Matrix instance to clone
     * @param {number} [nColumns] - Number of columns of the new matrix
     */
    class Matrix extends superCtor {
        static get [Symbol.species]() {
            return this;
        }

        /**
         * Constructs a Matrix with the chosen dimensions from a 1D array
         * @param {number} newRows - Number of rows
         * @param {number} newColumns - Number of columns
         * @param {Array} newData - A 1D array containing data for the matrix
         * @return {Matrix} - The new matrix
         */
        static from1DArray(newRows, newColumns, newData) {
            let length = newRows * newColumns;
            if (length !== newData.length) {
                throw new RangeError('Data length does not match given dimensions');
            }
            let newMatrix = new this(newRows, newColumns);
            for (let row = 0; row < newRows; row++) {
                for (let column = 0; column < newColumns; column++) {
                    newMatrix.set(row, column, newData[row * newColumns + column]);
                }
            }
            return newMatrix;
        }

        /**
         * Creates a row vector, a matrix with only one row.
         * @param {Array} newData - A 1D array containing data for the vector
         * @return {Matrix} - The new matrix
         */
        static rowVector(newData) {
            let vector = new this(1, newData.length);
            for (let i = 0; i < newData.length; i++) {
                vector.set(0, i, newData[i]);
            }
            return vector;
        }

        /**
         * Creates a column vector, a matrix with only one column.
         * @param {Array} newData - A 1D array containing data for the vector
         * @return {Matrix} - The new matrix
         */
        static columnVector(newData) {
            let vector = new this(newData.length, 1);
            for (let i = 0; i < newData.length; i++) {
                vector.set(i, 0, newData[i]);
            }
            return vector;
        }

        /**
         * Creates an empty matrix with the given dimensions. Values will be undefined. Same as using new Matrix(rows, columns).
         * @param {number} rows - Number of rows
         * @param {number} columns - Number of columns
         * @return {Matrix} - The new matrix
         */
        static empty(rows, columns) {
            return new this(rows, columns);
        }

        /**
         * Creates a matrix with the given dimensions. Values will be set to zero.
         * @param {number} rows - Number of rows
         * @param {number} columns - Number of columns
         * @return {Matrix} - The new matrix
         */
        static zeros(rows, columns) {
            return this.empty(rows, columns).fill(0);
        }

        /**
         * Creates a matrix with the given dimensions. Values will be set to one.
         * @param {number} rows - Number of rows
         * @param {number} columns - Number of columns
         * @return {Matrix} - The new matrix
         */
        static ones(rows, columns) {
            return this.empty(rows, columns).fill(1);
        }

        /**
         * Creates a matrix with the given dimensions. Values will be randomly set.
         * @param {number} rows - Number of rows
         * @param {number} columns - Number of columns
         * @param {function} [rng=Math.random] - Random number generator
         * @return {Matrix} The new matrix
         */
        static rand(rows, columns, rng) {
            if (rng === undefined) rng = Math.random;
            let matrix = this.empty(rows, columns);
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < columns; j++) {
                    matrix.set(i, j, rng());
                }
            }
            return matrix;
        }

        /**
         * Creates a matrix with the given dimensions. Values will be random integers.
         * @param {number} rows - Number of rows
         * @param {number} columns - Number of columns
         * @param {number} [minValue=0] - Minimal value
         * @param {number} [maxValue=1000] - Maximum value
         * @param {number} [decimal=5] - Decimal places after the point
         * @param {function} [rng=Math.random] - Random number generator
         * @return {Matrix} The new matrix
         */
        static randInt(rows, columns, minValue, maxValue, decimal, rng) {
            if (maxValue === undefined) maxValue = 1000;
            if (minValue === undefined) minValue = 0;
            if (decimal === undefined)  decimal = 5;
            if (rng === undefined) rng = Math.random;

            let matrix = this.empty(rows, columns);

            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < columns; j++) {
                    let value = Math.floor(rng() * (maxValue - minValue) + minValue).toFixed(decimal);
                    matrix.set(i, j, value);
                }
            }
            return matrix;
        }

        /**
         * Creates an identity matrix with the given dimension. Values of the diagonal will be 1 and others will be 0.
         * @param {number} rows - Number of rows
         * @param {number} [columns=rows] - Number of columns
         * @param {number} [value=1] - Value to fill the diagonal with
         * @return {Matrix} - The new identity matrix
         */
        static eye(rows, columns, value) {
            if (columns === undefined) columns = rows;
            if (value === undefined) value = 1;
            let min = Math.min(rows, columns);
            let matrix = this.zeros(rows, columns);
            for (let i = 0; i < min; i++) {
                matrix.set(i, i, value);
            }
            return matrix;
        }

        /**
         * Creates a diagonal matrix based on the given array.
         * @param {Array} data - Array containing the data for the diagonal
         * @param {number} [rows] - Number of rows (Default: data.length)
         * @param {number} [columns] - Number of columns (Default: rows)
         * @return {Matrix} - The new diagonal matrix
         */
        static diag(data, rows, columns) {
            let l = data.length;
            if (rows === undefined) rows = l;
            if (columns === undefined) columns = rows;
            let min = Math.min(l, rows, columns);
            let matrix = this.zeros(rows, columns);
            for (let i = 0; i < min; i++) {
                matrix.set(i, i, data[i]);
            }
            return matrix;
        }

        /**
         * Returns a matrix whose elements are the minimum between matrix1 and matrix2
         * @param {Matrix} matrix1
         * @param {Matrix} matrix2
         * @return {Matrix}
         */
        static min(matrix1, matrix2) {
            matrix1 = this.checkMatrix(matrix1);
            matrix2 = this.checkMatrix(matrix2);
            let rows = matrix1.rows;
            let columns = matrix1.columns;
            let result = new this(rows, columns);
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < columns; j++) {
                    result.set(i, j, Math.min(matrix1.get(i, j), matrix2.get(i, j)));
                }
            }
            return result;
        }

        /**
         * Returns a matrix whose elements are the maximum between matrix1 and matrix2
         * @param {Matrix} matrix1
         * @param {Matrix} matrix2
         * @return {Matrix}
         */
        static max(matrix1, matrix2) {
            matrix1 = this.checkMatrix(matrix1);
            matrix2 = this.checkMatrix(matrix2);
            let rows = matrix1.rows;
            let columns = matrix1.columns;
            let result = new this(rows, columns);
            for (let i = 0; i < rows; i++) {
                for (let j = 0; j < columns; j++) {
                    result.set(i, j, Math.max(matrix1.get(i, j), matrix2.get(i, j)));
                }
            }
            return result;
        }

        /**
         * Check that the provided value is a Matrix and tries to instantiate one if not
         * @param {*} value - The value to check
         * @return {Matrix}
         */
        static checkMatrix(value) {
            return Matrix.isMatrix(value) ? value : new this(value);
        }

        /**
         * Returns true if the argument is a Matrix, false otherwise
         * @param {*} value - The value to check
         * @return {boolean}
         */
        static isMatrix(value) {
            return (value != null) && (value.klass === 'Matrix');
        }

        /**
         * @prop {number} size - The number of elements in the matrix.
         */
        get size() {
            return this.rows * this.columns;
        }

        /**
         * Applies a callback for each element of the matrix. The function is called in the matrix (this) context.
         * @param {function} callback - Function that will be called with two parameters : i (row) and j (column)
         * @return {Matrix} this
         */
        apply(callback) {
            if (typeof callback !== 'function') {
                throw new TypeError('callback must be a function');
            }
            let ii = this.rows;
            let jj = this.columns;
            for (let i = 0; i < ii; i++) {
                for (let j = 0; j < jj; j++) {
                    callback.call(this, i, j);
                }
            }
            return this;
        }

        /**
         * Returns a new 1D array filled row by row with the matrix values
         * @return {Array}
         */
        to1DArray() {
            let array = new Array(this.size);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    array[i * this.columns + j] = this.get(i, j);
                }
            }
            return array;
        }

        /**
         * Returns a 2D array containing a copy of the data
         * @return {Array}
         */
        to2DArray() {
            let copy = new Array(this.rows);
            for (let i = 0; i < this.rows; i++) {
                copy[i] = new Array(this.columns);
                for (let j = 0; j < this.columns; j++) {
                    copy[i][j] = this.get(i, j);
                }
            }
            return copy;
        }

        /**
         * @return {boolean} true if the matrix has one row
         */
        isRowVector() {
            return this.rows === 1;
        }

        /**
         * @return {boolean} true if the matrix has one column
         */
        isColumnVector() {
            return this.columns === 1;
        }

        /**
         * @return {boolean} true if the matrix has one row or one column
         */
        isVector() {
            return (this.rows === 1) || (this.columns === 1);
        }

        /**
         * @return {boolean} true if the matrix has the same number of rows and columns
         */
        isSquare() {
            return this.rows === this.columns;
        }

        /**
         * @return {boolean} true if the matrix is square and has the same values on both sides of the diagonal
         */
        isSymmetric() {
            if (this.isSquare()) {
                for (let i = 0; i < this.rows; i++) {
                    for (let j = 0; j <= i; j++) {
                        if (this.get(i, j) !== this.get(j, i)) {
                            return false;
                        }
                    }
                }
                return true;
            }
            return false;
        }

        /**
         * Sets a given element of the matrix. mat.set(3,4,1) is equivalent to mat[3][4]=1
         * @abstract
         * @param {number} rowIndex - Index of the row
         * @param {number} columnIndex - Index of the column
         * @param {number} value - The new value for the element
         * @return {Matrix} this
         */
        set(rowIndex, columnIndex, value) { // eslint-disable-line no-unused-lets
            throw new Error('set method is unimplemented');
        }

        /**
         * Returns the given element of the matrix. mat.get(3,4) is equivalent to matrix[3][4]
         * @abstract
         * @param {number} rowIndex - Index of the row
         * @param {number} columnIndex - Index of the column
         * @return {number}
         */
        get(rowIndex, columnIndex) { // eslint-disable-line no-unused-lets
            throw new Error('get method is unimplemented');
        }

        /**
         * Creates a new matrix that is a repetition of the current matrix. New matrix has rowRep times the number of
         * rows of the matrix, and colRep times the number of columns of the matrix
         * @param {number} rowRep - Number of times the rows should be repeated
         * @param {number} colRep - Number of times the columns should be re
         * @return {Matrix}
         * @example
         * let matrix = new Matrix([[1,2]]);
         * matrix.repeat(2); // [[1,2],[1,2]]
         */
        repeat(rowRep, colRep) {
            rowRep = rowRep || 1;
            colRep = colRep || 1;
            let matrix = new this.constructor[Symbol.species](this.rows * rowRep, this.columns * colRep);
            for (let i = 0; i < rowRep; i++) {
                for (let j = 0; j < colRep; j++) {
                    matrix.setSubMatrix(this, this.rows * i, this.columns * j);
                }
            }
            return matrix;
        }

        /**
         * Fills the matrix with a given value. All elements will be set to this value.
         * @param {number} value - New value
         * @return {Matrix} this
         */
        fill(value) {
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, value);
                }
            }
            return this;
        }

        /**
         * Negates the matrix. All elements will be multiplied by (-1)
         * @return {Matrix} this
         */
        neg() {
            return this.mulS(-1);
        }

        /**
         * Returns a new array from the given row index
         * @param {number} index - Row index
         * @return {Array}
         */
        getRow(index) {
            checkRowIndex(this, index);
            let row = new Array(this.columns);
            for (let i = 0; i < this.columns; i++) {
                row[i] = this.get(index, i);
            }
            return row;
        }

        /**
         * Returns a new row vector from the given row index
         * @param {number} index - Row index
         * @return {Matrix}
         */
        getRowVector(index) {
            return this.constructor.rowVector(this.getRow(index));
        }

        /**
         * Sets a row at the given index
         * @param {number} index - Row index
         * @param {Array|Matrix} array - Array or vector
         * @return {Matrix} this
         */
        setRow(index, array) {
            checkRowIndex(this, index);
            array = checkRowVector(this, array);
            for (let i = 0; i < this.columns; i++) {
                this.set(index, i, array[i]);
            }
            return this;
        }

        /**
         * Swaps two rows
         * @param {number} row1 - First row index
         * @param {number} row2 - Second row index
         * @return {Matrix} this
         */
        swapRows(row1, row2) {
            checkRowIndex(this, row1);
            checkRowIndex(this, row2);
            for (let i = 0; i < this.columns; i++) {
                let temp = this.get(row1, i);
                this.set(row1, i, this.get(row2, i));
                this.set(row2, i, temp);
            }
            return this;
        }

        /**
         * Returns a new array from the given column index
         * @param {number} index - Column index
         * @return {Array}
         */
        getColumn(index) {
            checkColumnIndex(this, index);
            let column = new Array(this.rows);
            for (let i = 0; i < this.rows; i++) {
                column[i] = this.get(i, index);
            }
            return column;
        }

        /**
         * Returns a new column vector from the given column index
         * @param {number} index - Column index
         * @return {Matrix}
         */
        getColumnVector(index) {
            return this.constructor.columnVector(this.getColumn(index));
        }

        /**
         * Sets a column at the given index
         * @param {number} index - Column index
         * @param {Array|Matrix} array - Array or vector
         * @return {Matrix} this
         */
        setColumn(index, array) {
            checkColumnIndex(this, index);
            array = checkColumnVector(this, array);
            for (let i = 0; i < this.rows; i++) {
                this.set(i, index, array[i]);
            }
            return this;
        }

        /**
         * Swaps two columns
         * @param {number} column1 - First column index
         * @param {number} column2 - Second column index
         * @return {Matrix} this
         */
        swapColumns(column1, column2) {
            checkColumnIndex(this, column1);
            checkColumnIndex(this, column2);
            for (let i = 0; i < this.rows; i++) {
                let temp = this.get(i, column1);
                this.set(i, column1, this.get(i, column2));
                this.set(i, column2, temp);
            }
            return this;
        }

        /**
         * Adds the values of a vector to each row
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        addRowVector(vector) {
            vector = checkRowVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) + vector[j]);
                }
            }
            return this;
        }

        /**
         * Subtracts the values of a vector from each row
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        subRowVector(vector) {
            vector = checkRowVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) - vector[j]);
                }
            }
            return this;
        }

        /**
         * Multiplies the values of a vector with each row
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        mulRowVector(vector) {
            vector = checkRowVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) * vector[j]);
                }
            }
            return this;
        }

        /**
         * Divides the values of each row by those of a vector
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        divRowVector(vector) {
            vector = checkRowVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) / vector[j]);
                }
            }
            return this;
        }

        /**
         * Adds the values of a vector to each column
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        addColumnVector(vector) {
            vector = checkColumnVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) + vector[i]);
                }
            }
            return this;
        }

        /**
         * Subtracts the values of a vector from each column
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        subColumnVector(vector) {
            vector = checkColumnVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) - vector[i]);
                }
            }
            return this;
        }

        /**
         * Multiplies the values of a vector with each column
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        mulColumnVector(vector) {
            vector = checkColumnVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) * vector[i]);
                }
            }
            return this;
        }

        /**
         * Divides the values of each column by those of a vector
         * @param {Array|Matrix} vector - Array or vector
         * @return {Matrix} this
         */
        divColumnVector(vector) {
            vector = checkColumnVector(this, vector);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    this.set(i, j, this.get(i, j) / vector[i]);
                }
            }
            return this;
        }

        /**
         * Multiplies the values of a row with a scalar
         * @param {number} index - Row index
         * @param {number} value
         * @return {Matrix} this
         */
        mulRow(index, value) {
            checkRowIndex(this, index);
            for (let i = 0; i < this.columns; i++) {
                this.set(index, i, this.get(index, i) * value);
            }
            return this;
        }

        /**
         * Multiplies the values of a column with a scalar
         * @param {number} index - Column index
         * @param {number} value
         * @return {Matrix} this
         */
        mulColumn(index, value) {
            checkColumnIndex(this, index);
            for (let i = 0; i < this.rows; i++) {
                this.set(i, index, this.get(i, index) * value);
            }
            return this;
        }

        /**
         * Returns the maximum value of the matrix
         * @return {number}
         */
        max() {
            let v = this.get(0, 0);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    if (this.get(i, j) > v) {
                        v = this.get(i, j);
                    }
                }
            }
            return v;
        }

        /**
         * Returns the index of the maximum value
         * @return {Array}
         */
        maxIndex() {
            let v = this.get(0, 0);
            let idx = [0, 0];
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    if (this.get(i, j) > v) {
                        v = this.get(i, j);
                        idx[0] = i;
                        idx[1] = j;
                    }
                }
            }
            return idx;
        }

        /**
         * Returns the minimum value of the matrix
         * @return {number}
         */
        min() {
            let v = this.get(0, 0);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    if (this.get(i, j) < v) {
                        v = this.get(i, j);
                    }
                }
            }
            return v;
        }

        /**
         * Returns the index of the minimum value
         * @return {Array}
         */
        minIndex() {
            let v = this.get(0, 0);
            let idx = [0, 0];
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    if (this.get(i, j) < v) {
                        v = this.get(i, j);
                        idx[0] = i;
                        idx[1] = j;
                    }
                }
            }
            return idx;
        }

        /**
         * Returns the maximum value of one row
         * @param {number} row - Row index
         * @return {number}
         */
        maxRow(row) {
            checkRowIndex(this, row);
            let v = this.get(row, 0);
            for (let i = 1; i < this.columns; i++) {
                if (this.get(row, i) > v) {
                    v = this.get(row, i);
                }
            }
            return v;
        }

        /**
         * Returns the index of the maximum value of one row
         * @param {number} row - Row index
         * @return {Array}
         */
        maxRowIndex(row) {
            checkRowIndex(this, row);
            let v = this.get(row, 0);
            let idx = [row, 0];
            for (let i = 1; i < this.columns; i++) {
                if (this.get(row, i) > v) {
                    v = this.get(row, i);
                    idx[1] = i;
                }
            }
            return idx;
        }

        /**
         * Returns the minimum value of one row
         * @param {number} row - Row index
         * @return {number}
         */
        minRow(row) {
            checkRowIndex(this, row);
            let v = this.get(row, 0);
            for (let i = 1; i < this.columns; i++) {
                if (this.get(row, i) < v) {
                    v = this.get(row, i);
                }
            }
            return v;
        }

        /**
         * Returns the index of the maximum value of one row
         * @param {number} row - Row index
         * @return {Array}
         */
        minRowIndex(row) {
            checkRowIndex(this, row);
            let v = this.get(row, 0);
            let idx = [row, 0];
            for (let i = 1; i < this.columns; i++) {
                if (this.get(row, i) < v) {
                    v = this.get(row, i);
                    idx[1] = i;
                }
            }
            return idx;
        }

        /**
         * Returns the maximum value of one column
         * @param {number} column - Column index
         * @return {number}
         */
        maxColumn(column) {
            checkColumnIndex(this, column);
            let v = this.get(0, column);
            for (let i = 1; i < this.rows; i++) {
                if (this.get(i, column) > v) {
                    v = this.get(i, column);
                }
            }
            return v;
        }

        /**
         * Returns the index of the maximum value of one column
         * @param {number} column - Column index
         * @return {Array}
         */
        maxColumnIndex(column) {
            checkColumnIndex(this, column);
            let v = this.get(0, column);
            let idx = [0, column];
            for (let i = 1; i < this.rows; i++) {
                if (this.get(i, column) > v) {
                    v = this.get(i, column);
                    idx[0] = i;
                }
            }
            return idx;
        }

        /**
         * Returns the minimum value of one column
         * @param {number} column - Column index
         * @return {number}
         */
        minColumn(column) {
            checkColumnIndex(this, column);
            let v = this.get(0, column);
            for (let i = 1; i < this.rows; i++) {
                if (this.get(i, column) < v) {
                    v = this.get(i, column);
                }
            }
            return v;
        }

        /**
         * Returns the index of the minimum value of one column
         * @param {number} column - Column index
         * @return {Array}
         */
        minColumnIndex(column) {
            checkColumnIndex(this, column);
            let v = this.get(0, column);
            let idx = [0, column];
            for (let i = 1; i < this.rows; i++) {
                if (this.get(i, column) < v) {
                    v = this.get(i, column);
                    idx[0] = i;
                }
            }
            return idx;
        }

        /**
         * Returns an array containing the diagonal values of the matrix
         * @return {Array}
         */
        diag() {
            let min = Math.min(this.rows, this.columns);
            let diag = new Array(min);
            for (let i = 0; i < min; i++) {
                diag[i] = this.get(i, i);
            }
            return diag;
        }

        /**
         * Returns the sum by the argument given, if no argument given,
         * it returns the sum of all elements of the matrix.
         * @param {string} by - sum by 'row' or 'column'.
         * @return {Matrix|number}
         */
        sum(by) {
            switch (by) {
                case 'row':
                    return sumByRow(this);
                case 'column':
                    return sumByColumn(this);
                default:
                    return sumAll(this);
            }
        }

        /**
         * Returns the mean of all elements of the matrix
         * @return {number}
         */
        mean() {
            return this.sum() / this.size;
        }

        /**
         * Returns the product of all elements of the matrix
         * @return {number}
         */
        prod() {
            let prod = 1;
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    prod *= this.get(i, j);
                }
            }
            return prod;
        }

        /**
         * Returns the norm of a matrix.
         * @param {string} type - "frobenius" (default) or "max" return resp. the Frobenius norm and the max norm.
         * @return {number}
         */
        norm(type = 'frobenius') {
            let result = 0;
            if (type === 'max') {
                return this.max();
            } else if (type === 'frobenius') {
                for (let i = 0; i < this.rows; i++) {
                    for (let j = 0; j < this.columns; j++) {
                        result = result + this.get(i, j) * this.get(i, j);
                    }
                }
                return Math.sqrt(result);
            } else {
                throw new RangeError(`unknown norm type: ${type}`);
            }
        }

        /**
         * Computes the cumulative sum of the matrix elements (in place, row by row)
         * @return {Matrix} this
         */
        cumulativeSum() {
            let sum = 0;
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    sum += this.get(i, j);
                    this.set(i, j, sum);
                }
            }
            return this;
        }

        /**
         * Computes the dot (scalar) product between the matrix and another
         * @param {Matrix} vector2 vector
         * @return {number}
         */
        dot(vector2) {
            if (Matrix.isMatrix(vector2)) vector2 = vector2.to1DArray();
            let vector1 = this.to1DArray();
            if (vector1.length !== vector2.length) {
                throw new RangeError('vectors do not have the same size');
            }
            let dot = 0;
            for (let i = 0; i < vector1.length; i++) {
                dot += vector1[i] * vector2[i];
            }
            return dot;
        }

        /**
         * Returns the matrix product between this and other
         * @param {Matrix} other
         * @return {Matrix}
         */
        mmul(other) {
            other = this.constructor.checkMatrix(other);
            if (this.columns !== other.rows) {
                // eslint-disable-next-line no-console
                console.warn('Number of columns of left matrix are not equal to number of rows of right matrix.');
            }

            let m = this.rows;
            let n = this.columns;
            let p = other.columns;

            let result = new this.constructor[Symbol.species](m, p);

            let Bcolj = new Array(n);
            for (let j = 0; j < p; j++) {
                for (let k = 0; k < n; k++) {
                    Bcolj[k] = other.get(k, j);
                }

                for (let i = 0; i < m; i++) {
                    let s = 0;
                    for (k = 0; k < n; k++) {
                        s += this.get(i, k) * Bcolj[k];
                    }

                    result.set(i, j, s);
                }
            }
            return result;
        }

        strassen2x2(other) {
            let result = new this.constructor[Symbol.species](2, 2);
            const a11 = this.get(0, 0);
            const b11 = other.get(0, 0);
            const a12 = this.get(0, 1);
            const b12 = other.get(0, 1);
            const a21 = this.get(1, 0);
            const b21 = other.get(1, 0);
            const a22 = this.get(1, 1);
            const b22 = other.get(1, 1);

            // Compute intermediate values.
            const m1 = (a11 + a22) * (b11 + b22);
            const m2 = (a21 + a22) * b11;
            const m3 = a11 * (b12 - b22);
            const m4 = a22 * (b21 - b11);
            const m5 = (a11 + a12) * b22;
            const m6 = (a21 - a11) * (b11 + b12);
            const m7 = (a12 - a22) * (b21 + b22);

            // Combine intermediate values into the output.
            const c00 = m1 + m4 - m5 + m7;
            const c01 = m3 + m5;
            const c10 = m2 + m4;
            const c11 = m1 - m2 + m3 + m6;

            result.set(0, 0, c00);
            result.set(0, 1, c01);
            result.set(1, 0, c10);
            result.set(1, 1, c11);
            return result;
        }

        strassen3x3(other) {
            let result = new this.constructor[Symbol.species](3, 3);

            const a00 = this.get(0, 0);
            const a01 = this.get(0, 1);
            const a02 = this.get(0, 2);
            const a10 = this.get(1, 0);
            const a11 = this.get(1, 1);
            const a12 = this.get(1, 2);
            const a20 = this.get(2, 0);
            const a21 = this.get(2, 1);
            const a22 = this.get(2, 2);

            const b00 = other.get(0, 0);
            const b01 = other.get(0, 1);
            const b02 = other.get(0, 2);
            const b10 = other.get(1, 0);
            const b11 = other.get(1, 1);
            const b12 = other.get(1, 2);
            const b20 = other.get(2, 0);
            const b21 = other.get(2, 1);
            const b22 = other.get(2, 2);

            const m1 = (a00 + a01 + a02 - a10 - a11 - a21 - a22) * b11;
            const m2 = (a00 - a10) * (-b01 + b11);
            const m3 = a11 * (-b00 + b01 + b10 - b11 - b12 - b20 + b22);
            const m4 = (-a00 + a10 + a11) * (b00 - b01 + b11);
            const m5 = (a10 + a11) * (-b00 + b01);
            const m6 = a00 * b00;
            const m7 = (-a00 + a20 + a21) * (b00 - b02 + b12);
            const m8 = (-a00 + a20) * (b02 - b12);
            const m9 = (a20 + a21) * (-b00 + b02);
            const m10 = (a00 + a01 + a02 - a11 - a12 - a20 - a21) * b12;
            const m11 = a21 * (-b00 + b02 + b10 - b11 - b12 - b20 + b21);
            const m12 = (-a02 + a21 + a22) * (b11 + b20 - b21);
            const m13 = (a02 - a22) * (b11 - b21);
            const m14 = a02 * b20;
            const m15 = (a21 + a22) * (-b20 + b21);
            const m16 = (-a02 + a11 + a12) * (b12 + b20 - b22);
            const m17 = (a02 - a12) * (b12 - b22);
            const m18 = (a11 + a12) * (-b20 + b22);
            const m19 = a01 * b10;
            const m20 = a12 * b21;
            const m21 = a10 * b02;
            const m22 = a20 * b01;
            const m23 = a22 * b22;

            const c00 = m6 + m14 + m19;
            const c01 = m1 + m4 + m5 + m6 + m12 + m14 + m15;
            const c02 = m6 + m7 + m9 + m10 + m14 + m16 + m18;
            const c10 = m2 + m3 + m4 + m6 + m14 + m16 + m17;
            const c11 = m2 + m4 + m5 + m6 + m20;
            const c12 = m14 + m16 + m17 + m18 + m21;
            const c20 = m6 + m7 + m8 + m11 + m12 + m13 + m14;
            const c21 = m12 + m13 + m14 + m15 + m22;
            const c22 = m6 + m7 + m8 + m9 + m23;

            result.set(0, 0, c00);
            result.set(0, 1, c01);
            result.set(0, 2, c02);
            result.set(1, 0, c10);
            result.set(1, 1, c11);
            result.set(1, 2, c12);
            result.set(2, 0, c20);
            result.set(2, 1, c21);
            result.set(2, 2, c22);
            return result;
        }

        /**
         * Returns the matrix product between x and y. More efficient than mmul(other) only when we multiply squared matrix and when the size of the matrix is > 1000.
         * @param {Matrix} y
         * @return {Matrix}
         */
        mmulStrassen(y) {
            let x = this.clone();
            let r1 = x.rows;
            let c1 = x.columns;
            let r2 = y.rows;
            let c2 = y.columns;
            if (c1 !== r2) {
                // eslint-disable-next-line no-console
                console.warn(`Multiplying ${r1} x ${c1} and ${r2} x ${c2} matrix: dimensions do not match.`);
            }

            // Put a matrix into the top left of a matrix of zeros.
            // `rows` and `cols` are the dimensions of the output matrix.
            function embed(mat, rows, cols) {
                let r = mat.rows;
                let c = mat.columns;
                if ((r === rows) && (c === cols)) {
                    return mat;
                } else {
                    let resultat = Matrix.zeros(rows, cols);
                    resultat = resultat.setSubMatrix(mat, 0, 0);
                    return resultat;
                }
            }


            // Make sure both matrices are the same size.
            // This is exclusively for simplicity:
            // this algorithm can be implemented with matrices of different sizes.

            let r = Math.max(r1, r2);
            let c = Math.max(c1, c2);
            x = embed(x, r, c);
            y = embed(y, r, c);

            // Our recursive multiplication function.
            function blockMult(a, b, rows, cols) {
                // For small matrices, resort to naive multiplication.
                if (rows <= 512 || cols <= 512) {
                    return a.mmul(b); // a is equivalent to this
                }

                // Apply dynamic padding.
                if ((rows % 2 === 1) && (cols % 2 === 1)) {
                    a = embed(a, rows + 1, cols + 1);
                    b = embed(b, rows + 1, cols + 1);
                } else if (rows % 2 === 1) {
                    a = embed(a, rows + 1, cols);
                    b = embed(b, rows + 1, cols);
                } else if (cols % 2 === 1) {
                    a = embed(a, rows, cols + 1);
                    b = embed(b, rows, cols + 1);
                }

                let halfRows = parseInt(a.rows / 2);
                let halfCols = parseInt(a.columns / 2);
                // Subdivide input matrices.
                let a11 = a.subMatrix(0, halfRows - 1, 0, halfCols - 1);
                let b11 = b.subMatrix(0, halfRows - 1, 0, halfCols - 1);

                let a12 = a.subMatrix(0, halfRows - 1, halfCols, a.columns - 1);
                let b12 = b.subMatrix(0, halfRows - 1, halfCols, b.columns - 1);

                let a21 = a.subMatrix(halfRows, a.rows - 1, 0, halfCols - 1);
                let b21 = b.subMatrix(halfRows, b.rows - 1, 0, halfCols - 1);

                let a22 = a.subMatrix(halfRows, a.rows - 1, halfCols, a.columns - 1);
                let b22 = b.subMatrix(halfRows, b.rows - 1, halfCols, b.columns - 1);

                // Compute intermediate values.
                let m1 = blockMult(Matrix.add(a11, a22), Matrix.add(b11, b22), halfRows, halfCols);
                let m2 = blockMult(Matrix.add(a21, a22), b11, halfRows, halfCols);
                let m3 = blockMult(a11, Matrix.sub(b12, b22), halfRows, halfCols);
                let m4 = blockMult(a22, Matrix.sub(b21, b11), halfRows, halfCols);
                let m5 = blockMult(Matrix.add(a11, a12), b22, halfRows, halfCols);
                let m6 = blockMult(Matrix.sub(a21, a11), Matrix.add(b11, b12), halfRows, halfCols);
                let m7 = blockMult(Matrix.sub(a12, a22), Matrix.add(b21, b22), halfRows, halfCols);

                // Combine intermediate values into the output.
                let c11 = Matrix.add(m1, m4);
                c11.sub(m5);
                c11.add(m7);
                let c12 = Matrix.add(m3, m5);
                let c21 = Matrix.add(m2, m4);
                let c22 = Matrix.sub(m1, m2);
                c22.add(m3);
                c22.add(m6);

                //Crop output to the desired size (undo dynamic padding).
                let resultat = Matrix.zeros(2 * c11.rows, 2 * c11.columns);
                resultat = resultat.setSubMatrix(c11, 0, 0);
                resultat = resultat.setSubMatrix(c12, c11.rows, 0);
                resultat = resultat.setSubMatrix(c21, 0, c11.columns);
                resultat = resultat.setSubMatrix(c22, c11.rows, c11.columns);
                return resultat.subMatrix(0, rows - 1, 0, cols - 1);
            }
            return blockMult(x, y, r, c);
        }

        /**
         * Returns a row-by-row scaled matrix
         * @param {number} [min=0] - Minimum scaled value
         * @param {number} [max=1] - Maximum scaled value
         * @return {Matrix} - The scaled matrix
         */
        scaleRows(min, max) {
            min = min === undefined ? 0 : min;
            max = max === undefined ? 1 : max;
            if (min >= max) {
                throw new RangeError('min should be strictly smaller than max');
            }
            let newMatrix = this.constructor.empty(this.rows, this.columns);
            for (let i = 0; i < this.rows; i++) {
                let scaled = rescale(this.getRow(i), {min, max});
                newMatrix.setRow(i, scaled);
            }
            return newMatrix;
        }

        /**
         * Returns a new column-by-column scaled matrix
         * @param {number} [min=0] - Minimum scaled value
         * @param {number} [max=1] - Maximum scaled value
         * @return {Matrix} - The new scaled matrix
         * @example
         * let matrix = new Matrix([[1,2],[-1,0]]);
         * let scaledMatrix = matrix.scaleColumns(); // [[1,1],[0,0]]
         */
        scaleColumns(min, max) {
            min = min === undefined ? 0 : min;
            max = max === undefined ? 1 : max;
            if (min >= max) {
                throw new RangeError('min should be strictly smaller than max');
            }
            let newMatrix = this.constructor.empty(this.rows, this.columns);
            for (let i = 0; i < this.columns; i++) {
                let scaled = rescale(this.getColumn(i), {
                    min: min,
                    max: max
                });
                newMatrix.setColumn(i, scaled);
            }
            return newMatrix;
        }


        /**
         * Returns the Kronecker product (also known as tensor product) between this and other
         * See https://en.wikipedia.org/wiki/Kronecker_product
         * @param {Matrix} other
         * @return {Matrix}
         */
        kroneckerProduct(other) {
            other = this.constructor.checkMatrix(other);

            let m = this.rows;
            let n = this.columns;
            let p = other.rows;
            let q = other.columns;

            let result = new this.constructor[Symbol.species](m * p, n * q);
            for (let i = 0; i < m; i++) {
                for (let j = 0; j < n; j++) {
                    for (let k = 0; k < p; k++) {
                        for (let l = 0; l < q; l++) {
                            result[p * i + k][q * j + l] = this.get(i, j) * other.get(k, l);
                        }
                    }
                }
            }
            return result;
        }

        /**
         * Transposes the matrix and returns a new one containing the result
         * @return {Matrix}
         */
        transpose() {
            let result = new this.constructor[Symbol.species](this.columns, this.rows);
            for (let i = 0; i < this.rows; i++) {
                for (let j = 0; j < this.columns; j++) {
                    result.set(j, i, this.get(i, j));
                }
            }
            return result;
        }

        /**
         * Sorts the rows (in place)
         * @param {function} compareFunction - usual Array.prototype.sort comparison function
         * @return {Matrix} this
         */
        sortRows(compareFunction) {
            if (compareFunction === undefined) compareFunction = compareNumbers;
            for (let i = 0; i < this.rows; i++) {
                this.setRow(i, this.getRow(i).sort(compareFunction));
            }
            return this;
        }

        /**
         * Sorts the columns (in place)
         * @param {function} compareFunction - usual Array.prototype.sort comparison function
         * @return {Matrix} this
         */
        sortColumns(compareFunction) {
            if (compareFunction === undefined) compareFunction = compareNumbers;
            for (let i = 0; i < this.columns; i++) {
                this.setColumn(i, this.getColumn(i).sort(compareFunction));
            }
            return this;
        }

        /**
         * Returns a subset of the matrix
         * @param {number} startRow - First row index
         * @param {number} endRow - Last row index
         * @param {number} startColumn - First column index
         * @param {number} endColumn - Last column index
         * @return {Matrix}
         */
        subMatrix(startRow, endRow, startColumn, endColumn) {
            checkRange(this, startRow, endRow, startColumn, endColumn);
            let newMatrix = new this.constructor[Symbol.species](endRow - startRow + 1, endColumn - startColumn + 1);
            for (let i = startRow; i <= endRow; i++) {
                for (let j = startColumn; j <= endColumn; j++) {
                    newMatrix[i - startRow][j - startColumn] = this.get(i, j);
                }
            }
            return newMatrix;
        }

        /**
         * Returns a subset of the matrix based on an array of row indices
         * @param {Array} indices - Array containing the row indices
         * @param {number} [startColumn = 0] - First column index
         * @param {number} [endColumn = this.columns-1] - Last column index
         * @return {Matrix}
         */
        subMatrixRow(indices, startColumn, endColumn) {
            if (startColumn === undefined) startColumn = 0;
            if (endColumn === undefined) endColumn = this.columns - 1;
            if ((startColumn > endColumn) || (startColumn < 0) || (startColumn >= this.columns) || (endColumn < 0) || (endColumn >= this.columns)) {
                throw new RangeError('Argument out of range');
            }

            let newMatrix = new this.constructor[Symbol.species](indices.length, endColumn - startColumn + 1);
            for (let i = 0; i < indices.length; i++) {
                for (let j = startColumn; j <= endColumn; j++) {
                    if (indices[i] < 0 || indices[i] >= this.rows) {
                        throw new RangeError('Row index out of range: ' + indices[i]);
                    }
                    newMatrix.set(i, j - startColumn, this.get(indices[i], j));
                }
            }
            return newMatrix;
        }

        /**
         * Returns a subset of the matrix based on an array of column indices
         * @param {Array} indices - Array containing the column indices
         * @param {number} [startRow = 0] - First row index
         * @param {number} [endRow = this.rows-1] - Last row index
         * @return {Matrix}
         */
        subMatrixColumn(indices, startRow, endRow) {
            if (startRow === undefined) startRow = 0;
            if (endRow === undefined) endRow = this.rows - 1;
            if ((startRow > endRow) || (startRow < 0) || (startRow >= this.rows) || (endRow < 0) || (endRow >= this.rows)) {
                throw new RangeError('Argument out of range');
            }

            let newMatrix = new this.constructor[Symbol.species](endRow - startRow + 1, indices.length);
            for (let i = 0; i < indices.length; i++) {
                for (let j = startRow; j <= endRow; j++) {
                    if (indices[i] < 0 || indices[i] >= this.columns) {
                        throw new RangeError('Column index out of range: ' + indices[i]);
                    }
                    newMatrix.set(j - startRow, i, this.get(j, indices[i]));
                }
            }
            return newMatrix;
        }

        /**
         * Set a part of the matrix to the given sub-matrix
         * @param {Matrix|Array< Array >} matrix - The source matrix from which to extract values.
         * @param {number} startRow - The index of the first row to set
         * @param {number} startColumn - The index of the first column to set
         * @return {Matrix}
         */
        setSubMatrix(matrix, startRow, startColumn) {
            matrix = this.constructor.checkMatrix(matrix);
            let endRow = startRow + matrix.rows - 1;
            let endColumn = startColumn + matrix.columns - 1;
            checkRange(this, startRow, endRow, startColumn, endColumn);
            for (let i = 0; i < matrix.rows; i++) {
                for (let j = 0; j < matrix.columns; j++) {
                    this[startRow + i][startColumn + j] = matrix.get(i, j);
                }
            }
            return this;
        }

        /**
         * Return a new matrix based on a selection of rows and columns
         * @param {Array<number>} rowIndices - The row indices to select. Order matters and an index can be more than once.
         * @param {Array<number>} columnIndices - The column indices to select. Order matters and an index can be use more than once.
         * @return {Matrix} The new matrix
         */
        selection(rowIndices, columnIndices) {
            let indices = checkIndices(this, rowIndices, columnIndices);
            let newMatrix = new this.constructor[Symbol.species](rowIndices.length, columnIndices.length);
            for (let i = 0; i < indices.row.length; i++) {
                let rowIndex = indices.row[i];
                for (let j = 0; j < indices.column.length; j++) {
                    let columnIndex = indices.column[j];
                    newMatrix[i][j] = this.get(rowIndex, columnIndex);
                }
            }
            return newMatrix;
        }

        /**
         * Returns the trace of the matrix (sum of the diagonal elements)
         * @return {number}
         */
        trace() {
            let min = Math.min(this.rows, this.columns);
            let trace = 0;
            for (let i = 0; i < min; i++) {
                trace += this.get(i, i);
            }
            return trace;
        }

        /*
         Matrix views
         */

        /**
         * Returns a view of the transposition of the matrix
         * @return {MatrixTransposeView}
         */
        transposeView() {
            return new MatrixTransposeView(this);
        }

        /**
         * Returns a view of the row vector with the given index
         * @param {number} row - row index of the vector
         * @return {MatrixRowView}
         */
        rowView(row) {
            checkRowIndex(this, row);
            return new MatrixRowView(this, row);
        }

        /**
         * Returns a view of the column vector with the given index
         * @param {number} column - column index of the vector
         * @return {MatrixColumnView}
         */
        columnView(column) {
            checkColumnIndex(this, column);
            return new MatrixColumnView(this, column);
        }

        /**
         * Returns a view of the matrix flipped in the row axis
         * @return {MatrixFlipRowView}
         */
        flipRowView() {
            return new MatrixFlipRowView(this);
        }

        /**
         * Returns a view of the matrix flipped in the column axis
         * @return {MatrixFlipColumnView}
         */
        flipColumnView() {
            return new MatrixFlipColumnView(this);
        }

        /**
         * Returns a view of a submatrix giving the index boundaries
         * @param {number} startRow - first row index of the submatrix
         * @param {number} endRow - last row index of the submatrix
         * @param {number} startColumn - first column index of the submatrix
         * @param {number} endColumn - last column index of the submatrix
         * @return {MatrixSubView}
         */
        subMatrixView(startRow, endRow, startColumn, endColumn) {
            return new MatrixSubView(this, startRow, endRow, startColumn, endColumn);
        }

        /**
         * Returns a view of the cross of the row indices and the column indices
         * @example
         * // resulting vector is [[2], [2]]
         * let matrix = new Matrix([[1,2,3], [4,5,6]]).selectionView([0, 0], [1])
         * @param {Array<number>} rowIndices
         * @param {Array<number>} columnIndices
         * @return {MatrixSelectionView}
         */
        selectionView(rowIndices, columnIndices) {
            return new MatrixSelectionView(this, rowIndices, columnIndices);
        }

        /**
         * Returns a view of the row indices
         * @example
         * // resulting vector is [[1,2,3], [1,2,3]]
         * let matrix = new Matrix([[1,2,3], [4,5,6]]).rowSelectionView([0, 0])
         * @param {Array<number>} rowIndices
         * @return {MatrixRowSelectionView}
         */
        rowSelectionView(rowIndices) {
            return new MatrixRowSelectionView(this, rowIndices);
        }

        /**
         * Returns a view of the column indices
         * @example
         * // resulting vector is [[2, 2], [5, 5]]
         * let matrix = new Matrix([[1,2,3], [4,5,6]]).columnSelectionView([1, 1])
         * @param {Array<number>} columnIndices
         * @return {MatrixColumnSelectionView}
         */
        columnSelectionView(columnIndices) {
            return new MatrixColumnSelectionView(this, columnIndices);
        }


        /**
        * Calculates and returns the determinant of a matrix as a Number
        * @example
        *   new Matrix([[1,2,3], [4,5,6]]).det()
        * @return {number}
        */
        det() {
            if (this.isSquare()) {
                let a, b, c, d;
                if (this.columns === 2) {
                    // 2 x 2 matrix
                    a = this.get(0, 0);
                    b = this.get(0, 1);
                    c = this.get(1, 0);
                    d = this.get(1, 1);

                    return a * d - (b * c);
                } else if (this.columns === 3) {
                    // 3 x 3 matrix
                    let subMatrix0, subMatrix1, subMatrix2;
                    subMatrix0 = this.selectionView([1, 2], [1, 2]);
                    subMatrix1 = this.selectionView([1, 2], [0, 2]);
                    subMatrix2 = this.selectionView([1, 2], [0, 1]);
                    a = this.get(0, 0);
                    b = this.get(0, 1);
                    c = this.get(0, 2);

                    return a * subMatrix0.det() - b * subMatrix1.det() + c * subMatrix2.det();
                } else {
                    // general purpose determinant using the LU decomposition
                    return new LuDecomposition(this).determinant;
                }

            } else {
                throw Error('Determinant can only be calculated for a square matrix.');
            }
        }

        /**
         * Returns inverse of a matrix if it exists or the pseudoinverse
         * @param {number} threshold - threshold for taking inverse of singular values (default = 1e-15)
         * @return {Matrix} the (pseudo)inverted matrix.
         */
        pseudoInverse(threshold) {
            if (threshold === undefined) threshold = Number.EPSILON;
            let svdSolution = new SvDecomposition(this, {autoTranspose: true});

            let U = svdSolution.leftSingularVectors;
            let V = svdSolution.rightSingularVectors;
            let s = svdSolution.diagonal;

            for (let i = 0; i < s.length; i++) {
                if (Math.abs(s[i]) > threshold) {
                    s[i] = 1.0 / s[i];
                } else {
                    s[i] = 0.0;
                }
            }

            // convert list to diagonal
            s = this.constructor[Symbol.species].diag(s);
            return V.mmul(s.mmul(U.transposeView()));
        }

        /**
         * Creates an exact and independent copy of the matrix
         * @return {Matrix}
         */
        clone() {
            let newMatrix = new this.constructor[Symbol.species](this.rows, this.columns);
            for (let row = 0; row < this.rows; row++) {
                for (let column = 0; column < this.columns; column++) {
                    newMatrix.set(row, column, this.get(row, column));
                }
            }
            return newMatrix;
        }
    }

    Matrix.prototype.klass = 'Matrix';

    function compareNumbers(a, b) {
        return a - b;
    }

    /*
     Synonyms
     */

    Matrix.random = Matrix.rand;
    Matrix.diagonal = Matrix.diag;
    Matrix.prototype.diagonal = Matrix.prototype.diag;
    Matrix.identity = Matrix.eye;
    Matrix.prototype.negate = Matrix.prototype.neg;
    Matrix.prototype.tensorProduct = Matrix.prototype.kroneckerProduct;
    Matrix.prototype.determinant = Matrix.prototype.det;

    /*
     Add dynamically instance and static methods for mathematical operations
     */

    let inplaceOperator = `
(function %name%(value) {
    if (typeof value === 'number') return this.%name%S(value);
    return this.%name%M(value);
})
`;

    let inplaceOperatorScalar = `
(function %name%S(value) {
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, this.get(i, j) %op% value);
        }
    }
    return this;
})
`;

    let inplaceOperatorMatrix = `
(function %name%M(matrix) {
    matrix = this.constructor.checkMatrix(matrix);
    if (this.rows !== matrix.rows ||
        this.columns !== matrix.columns) {
        throw new RangeError('Matrices dimensions must be equal');
    }
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, this.get(i, j) %op% matrix.get(i, j));
        }
    }
    return this;
})
`;

    let staticOperator = `
(function %name%(matrix, value) {
    let newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%(value);
})
`;

    let inplaceMethod = `
(function %name%() {
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j)));
        }
    }
    return this;
})
`;

    let staticMethod = `
(function %name%(matrix) {
    let newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%();
})
`;

    let inplaceMethodWithArgs = `
(function %name%(%args%) {
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), %args%));
        }
    }
    return this;
})
`;

    let staticMethodWithArgs = `
(function %name%(matrix, %args%) {
    let newMatrix = new this[Symbol.species](matrix);
    return newMatrix.%name%(%args%);
})
`;


    let inplaceMethodWithOneArgScalar = `
(function %name%S(value) {
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), value));
        }
    }
    return this;
})
`;
    let inplaceMethodWithOneArgMatrix = `
(function %name%M(matrix) {
    matrix = this.constructor.checkMatrix(matrix);
    if (this.rows !== matrix.rows ||
        this.columns !== matrix.columns) {
        throw new RangeError('Matrices dimensions must be equal');
    }
    for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.columns; j++) {
            this.set(i, j, %method%(this.get(i, j), matrix.get(i, j)));
        }
    }
    return this;
})
`;

    let inplaceMethodWithOneArg = `
(function %name%(value) {
    if (typeof value === 'number') return this.%name%S(value);
    return this.%name%M(value);
})
`;

    let staticMethodWithOneArg = staticMethodWithArgs;

    let operators = [
        // Arithmetic operators
        ['+', 'add'],
        ['-', 'sub', 'subtract'],
        ['*', 'mul', 'multiply'],
        ['/', 'div', 'divide'],
        ['%', 'mod', 'modulus'],
        // Bitwise operators
        ['&', 'and'],
        ['|', 'or'],
        ['^', 'xor'],
        ['<<', 'leftShift'],
        ['>>', 'signPropagatingRightShift'],
        ['>>>', 'rightShift', 'zeroFillRightShift']
    ];

    let i;
    let eval2 = eval;
    for (let operator of operators) {
        let inplaceOp = eval2(fillTemplateFunction(inplaceOperator, {name: operator[1], op: operator[0]}));
        let inplaceOpS = eval2(fillTemplateFunction(inplaceOperatorScalar, {name: operator[1] + 'S', op: operator[0]}));
        let inplaceOpM = eval2(fillTemplateFunction(inplaceOperatorMatrix, {name: operator[1] + 'M', op: operator[0]}));
        let staticOp = eval2(fillTemplateFunction(staticOperator, {name: operator[1]}));
        for (i = 1; i < operator.length; i++) {
            Matrix.prototype[operator[i]] = inplaceOp;
            Matrix.prototype[operator[i] + 'S'] = inplaceOpS;
            Matrix.prototype[operator[i] + 'M'] = inplaceOpM;
            Matrix[operator[i]] = staticOp;
        }
    }

    let methods = [
        ['~', 'not']
    ];

    [
        'abs', 'acos', 'acosh', 'asin', 'asinh', 'atan', 'atanh', 'cbrt', 'ceil',
        'clz32', 'cos', 'cosh', 'exp', 'expm1', 'floor', 'fround', 'log', 'log1p',
        'log10', 'log2', 'round', 'sign', 'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'trunc'
    ].forEach(function (mathMethod) {
        methods.push(['Math.' + mathMethod, mathMethod]);
    });

    for (let method of methods) {
        let inplaceMeth = eval2(fillTemplateFunction(inplaceMethod, {name: method[1], method: method[0]}));
        let staticMeth = eval2(fillTemplateFunction(staticMethod, {name: method[1]}));
        for (i = 1; i < method.length; i++) {
            Matrix.prototype[method[i]] = inplaceMeth;
            Matrix[method[i]] = staticMeth;
        }
    }

    let methodsWithArgs = [
        ['Math.pow', 1, 'pow']
    ];

    for (let methodWithArg of methodsWithArgs) {
        let args = 'arg0';
        for (i = 1; i < methodWithArg[1]; i++) {
            args += `, arg${i}`;
        }
        if (methodWithArg[1] !== 1) {
            let inplaceMethWithArgs = eval2(fillTemplateFunction(inplaceMethodWithArgs, {
                name: methodWithArg[2],
                method: methodWithArg[0],
                args: args
            }));
            let staticMethWithArgs = eval2(fillTemplateFunction(staticMethodWithArgs, {name: methodWithArg[2], args: args}));
            for (i = 2; i < methodWithArg.length; i++) {
                Matrix.prototype[methodWithArg[i]] = inplaceMethWithArgs;
                Matrix[methodWithArg[i]] = staticMethWithArgs;
            }
        } else {
            let tmpllet = {
                name: methodWithArg[2],
                args: args,
                method: methodWithArg[0]
            };
            let inplaceMethod2 = eval2(fillTemplateFunction(inplaceMethodWithOneArg, tmpllet));
            let inplaceMethodS = eval2(fillTemplateFunction(inplaceMethodWithOneArgScalar, tmpllet));
            let inplaceMethodM = eval2(fillTemplateFunction(inplaceMethodWithOneArgMatrix, tmpllet));
            let staticMethod2 = eval2(fillTemplateFunction(staticMethodWithOneArg, tmpllet));
            for (i = 2; i < methodWithArg.length; i++) {
                Matrix.prototype[methodWithArg[i]] = inplaceMethod2;
                Matrix.prototype[methodWithArg[i] + 'M'] = inplaceMethodM;
                Matrix.prototype[methodWithArg[i] + 'S'] = inplaceMethodS;
                Matrix[methodWithArg[i]] = staticMethod2;
            }
        }
    }

    function fillTemplateFunction(template, values) {
        for (let value in values) {
            template = template.replace(new RegExp('%' + value + '%', 'g'), values[value]);
        }
        return template;
    }

    return Matrix;
}
