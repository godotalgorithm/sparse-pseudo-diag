// minimal C implementation of sparse pseudo-diagonalization as in the fast MOZYME solver of MOPAC
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// This implementation favors simplicity over performance, although it retains the essential structure
// needed for reasonable baseline performance. MOZYME improves performance by making use of atom-based
// block matrix structure, symmetry reduction of matrix storage (lower triangle), and geometric
// information for screening of LMO-LMO interactions. Empirical comparisons suggest that MOZYME also
// has additional mechanisms for pruning the list of occupied-virtual rotations which are necessary to
// maintain good sparsity in 3D protein structures that have not yet been clearly identified and
// replicated here. Bridging this gap is the subject of ongoing inquiry.

// The three numerical parameters below are set to reasonable, somewhat aggressive values that are
// roughly comparable to typical MOZYME calculations of large organic molecules.

// default cutoff value for removing matrix elements from the LMO matrix
#define LMO_CUTOFF 1e-7

// default cutoff value for removing matrix elements from the Fock matrix in the LMO basis
#define FOCK_CUTOFF 1e-7

// default value of convergence tolerance for Fock matrix elements coupling occupied & virtual LMOs
#define FOCK_TOLERANCE 1e-2

// buffered sparse vector
struct sparse_vec
{
    int size, num; // size of buffer & # of nonzero vector elements
    int *row; // row index of vector elements [size]
    double *val; // value of vector elements [size]
};

// square sparse matrix structure: similar to compressed-sparse-column format w/ data splitting over columns
struct sparse_mat
{
    int ncol; // number of columns
    struct sparse_vec *col; // sparse column vectors of the matrix [ncol]
};

// basic allocation/deallocation functions for sparse vectors & matrices
void malloc_vec(int size, struct sparse_vec *vec)
{
    vec->size = size;
    vec->row = (int*)malloc(sizeof(int)*size);
    vec->val = (double*)malloc(sizeof(double)*size);
}
void realloc_vec(int size, struct sparse_vec *vec)
{
    vec->size = size;
    vec->row = (int*)realloc(vec->row, sizeof(int)*size);
    vec->val = (double*)realloc(vec->val, sizeof(double)*size);
}
void free_vec(struct sparse_vec *vec)
{ if(vec->size > 0) { free(vec->row); free(vec->val); } }
void free_mat(struct sparse_mat *mat)
{ for(int i=0 ; i<mat->ncol ; i++) { if(mat->col[i].size > 0) { free_vec(mat->col+i); } } free(mat->col); }

// read sparse matrix from a text file (ncol preset in mat)
struct pair { int i; double d; };
int cmp2(const void *a, const void *b) { return ( ((struct pair*)a)->i - ((struct pair*)b)->i ); }
void read_mat(struct sparse_mat *mat, FILE *data)
{
    mat->col = (struct sparse_vec*)malloc(sizeof(struct sparse_vec)*mat->ncol);
    struct pair *buffer = (struct pair*)malloc(sizeof(struct pair)*mat->ncol);
    for(int i=0 ; i<mat->ncol ; i++)
    {
        int size;
        fscanf(data, "%d", &size);
        malloc_vec(size, mat->col+i);
        // read columns into a temporary buffer to use qsort to sort rows, which don't have guaranteed order
        for(int j=0 ; j<size ; j++)
        {
            fscanf(data, "%d %lf", &buffer[j].i, &buffer[j].d);
            buffer[j].i--; // switch from 1-based indexing to 0-based indexing
        }
        qsort(buffer, size, sizeof(struct pair), cmp2);
        for(int j=0 ; j<size ; j++)
        {
            mat->col[i].row[j] = buffer[j].i;
            mat->col[i].val[j] = buffer[j].d;
        }
        mat->col[i].num = mat->col[i].size;
    }
    free(buffer);
}

// build the transpose of a sparse matrix (ncol preset in out)
void transpose_mat(int nrow, struct sparse_mat *mat, struct sparse_mat *trans)
{
    trans->ncol = nrow;
    trans->col = (struct sparse_vec*)malloc(sizeof(struct sparse_vec)*trans->ncol);
    for(int i=0 ; i<trans->ncol ; i++)
    { trans->col[i].size = trans->col[i].num = 0; }
    // 1st pass: determine # of matrix elements per column
    for(int i=0 ; i<mat->ncol ; i++)
    for(int j=0 ; j<mat->col[i].num ; j++)
    { trans->col[mat->col[i].row[j]].size++; }
    for(int i=0 ; i<trans->ncol ; i++)
    { malloc_vec(trans->col[i].size, trans->col+i); }
    // 2nd pass: fill in matrix elements
    for(int i=0 ; i<mat->ncol ; i++)
    for(int j=0 ; j<mat->col[i].num ; j++)
    {
        int k = mat->col[i].row[j];
        trans->col[k].row[trans->col[k].num] = i;
        trans->col[k].val[trans->col[k].num++] = mat->col[i].val[j];
    }
}

// add two sparse vectors
void add_vec(double wt1, double wt2, struct sparse_vec *in1, struct sparse_vec *in2, struct sparse_vec *out)
{
    out->num = 0;
    int i1 = 0, i2 = 0;
    // merge vectors
    while(i1 < in1->num || i2 < in2->num)
    {
        if(out->num == out->size)
        { printf("ERROR: vector buffer overflow\n"); exit(1); }
        if(i1 != in1->num && (i2 == in2->num || in1->row[i1] < in2->row[i2]))
        {
            out->row[out->num] = in1->row[i1];
            out->val[out->num++] = wt1*in1->val[i1++];
        }
        else if(i1 == in1->num || in1->row[i1] > in2->row[i2])
        {
            out->row[out->num] = in2->row[i2];
            out->val[out->num++] = wt2*in2->val[i2++];
        }
        else
        {
            out->row[out->num] = in1->row[i1];
            out->val[out->num++] = wt1*in1->val[i1++] + wt2*in2->val[i2++];
        }
    }
}

// sparse matrix / sparse vector multiplication
void mat_vec(struct sparse_mat *mat, struct sparse_vec *in, struct sparse_vec *out, struct sparse_vec *work)
{
    // extra swap as needed to guarantee output ends in out
    if(in->num%2) { struct sparse_vec *swap = out; out = work; work = swap; }
    out->num = 0;
    for(int i=0 ; i<in->num ; i++)
    {
        struct sparse_vec *swap = out; out = work; work = swap;
        add_vec(1.0, in->val[i], work, mat->col+in->row[i], out);
    }
}

// copy a sparse vector & reallocate memory as needed
void copy_vec(struct sparse_vec *vec, struct sparse_vec *copy)
{
    if(copy->size < vec->num)
    {
        free_vec(copy);
        malloc_vec(vec->num, copy);
    }
    copy->num = vec->num;
    for(int i=0 ; i<copy->num ; i++)
    {
        copy->row[i] = vec->row[i];
        copy->val[i] = vec->val[i];
    }
}

// remove small elements from a sparse vector
void prune_vec(double cut, struct sparse_vec *vec)
{
    int num_cut = 0;
    for(int i=0 ; i<vec->num ; i++)
    {
        if(fabs(vec->val[i]) <= cut)
        { num_cut++; continue; }
        vec->row[i-num_cut] = vec->row[i];
        vec->val[i-num_cut] = vec->val[i];
    }
    vec->num -= num_cut;
}

// dot product of two sparse vectors
double dot_vec(struct sparse_vec *vec1, struct sparse_vec *vec2)
{
    double prod = 0.0;
    int i1 = 0, i2 = 0;
    while(i1 < vec1->num && i2 < vec2->num)
    {
        if(vec1->row[i1] == vec2->row[i2])
        { prod += vec1->val[i1++]*vec2->val[i2++]; }
        else if(vec1->row[i1] < vec2->row[i2])
        { i1++; }
        else
        { i2++; }
    }
    return prod;
}

// project Fock matrix into LMO basis & recompute occupied energies
void project_fock(double cut, struct sparse_mat *fock, struct sparse_mat *occ, struct sparse_mat *vir_trans, struct sparse_mat *fock_lmo)
{
    struct sparse_vec work, work2, work3;
    malloc_vec(fock->ncol, &work); malloc_vec(fock->ncol, &work2); malloc_vec(fock->ncol, &work3);
    fock_lmo->ncol = occ->ncol;
    fock_lmo->col = (struct sparse_vec*)malloc(sizeof(struct sparse_vec)*fock_lmo->ncol);
    for(int i=0 ; i<fock_lmo->ncol ; i++)
    { fock_lmo->col[i].size = fock_lmo->col[i].num = 0; }
    // form projected Fock matrix, one column at a time
    for(int i=0 ; i<fock_lmo->ncol ; i++)
    {
        mat_vec(fock, occ->col+i, &work, &work2);
        mat_vec(vir_trans, &work, &work2, &work3);
        prune_vec(cut, &work2);
        copy_vec(&work2, fock_lmo->col+i);
    }
    free_vec(&work); free_vec(&work2); free_vec(&work3);
}

// compute LMO energies
void project_energy(struct sparse_mat *fock, struct sparse_mat *lmo, double *energy)
{
    struct sparse_vec work, work2;
    malloc_vec(fock->ncol, &work); malloc_vec(fock->ncol, &work2);
    for(int i=0 ; i<lmo->ncol ; i++)
    {
        mat_vec(fock, lmo->col+i, &work, &work2);
        energy[i] = dot_vec(lmo->col+i, &work);
    }
    free_vec(&work); free_vec(&work2);
}

// scale a vector
void scale_vec(double wt, struct sparse_vec *vec)
{
    for(int i=0 ; i<vec->num ; i++)
    { vec->val[i] *= wt; }
}

// renormalize columns of a matrix
void renormalize(struct sparse_mat *mat)
{
    for(int i=0 ; i<mat->ncol ; i++)
    {
        double norm = 1.0/sqrt(dot_vec(mat->col+i, mat->col+i));
        scale_vec(norm, mat->col+i);
    }
}

// replace a row in a sparse matrix
int cmp(const void *a, const void *b) { return ( *(int*)a - *(int*)b ); }
void copy_row(int index, struct sparse_vec *old_vec, struct sparse_vec *new_vec, struct sparse_mat *mat)
{
    int iold=0, inew=0;
    while(iold < old_vec->num || inew < new_vec->num)
    {
        // insertion
        if(inew != new_vec->num && (iold == old_vec->num || old_vec->row[iold] > new_vec->row[inew]))
        {
            int i;
            struct sparse_vec *target = mat->col + new_vec->row[inew];
            if(target->num == target->size)
            { realloc_vec(++target->size, target); }
            for(i=target->num ; i>0 && target->row[i-1] > index ; i--)
            {
                target->row[i] = target->row[i-1];
                target->val[i] = target->val[i-1];
            }
            target->row[i] = index;
            target->val[i] = new_vec->val[inew];
            target->num++;
            inew++;
        }
        // removal
        else if(inew == new_vec->num || old_vec->row[iold] < new_vec->row[inew])
        {
            int offset = 0;
            struct sparse_vec *target = mat->col + old_vec->row[iold];
            for(int i=0 ; i<target->num-1 ; i++)
            {
                if(target->row[i] == index)
                { offset++; }
                target->row[i] = target->row[i+offset];
                target->val[i] = target->val[i+offset];
            }
            target->num--;
            iold++;
        }
        // replacement
        else
        {
            struct sparse_vec *target = mat->col + new_vec->row[inew];
            int offset = (int*)bsearch(&index, target->row, target->num, sizeof(int), cmp) - target->row;
            target->val[offset] = new_vec->val[inew];
            iold++; inew++;
        }
    }
}

// Jacobi rotation w/ sparsification applied to LMO & Fock matrices
void jacobi(double cut, double rot, int iocc, int ivir, struct sparse_mat *occ, struct sparse_mat *vir, 
            struct sparse_mat *occ_trans, struct sparse_mat *vir_trans, struct sparse_vec *work, struct sparse_vec *work2)
{
    // form updated LMO columns
    double c = 1.0/sqrt(1.0 + rot*rot);
    double s = c*rot;
    add_vec(c, -s, occ->col+iocc, vir->col+ivir, work);
    add_vec(s, c, occ->col+iocc, vir->col+ivir, work2);
    // remove small matrix elements
    prune_vec(cut, work);
    prune_vec(cut, work2);
    // renormalize
    double norm = 1.0/sqrt(dot_vec(work, work));
    scale_vec(norm, work);
    norm = 1.0/sqrt(dot_vec(work2, work2));
    scale_vec(norm, work2);
    // apply updates to transposes
    copy_row(iocc, occ->col+iocc, work, occ_trans);
    copy_row(ivir, vir->col+ivir, work2, vir_trans);
    // apply updates to matrices
    copy_vec(work, occ->col+iocc);
    copy_vec(work2, vir->col+ivir);
}

// pseudo-diagonalization solver
void pseudo_diag(double cut_lmo, double cut_fock, double tol, double *energy_occ, double *energy_vir, struct sparse_mat *fock,
                 struct sparse_mat *occ, struct sparse_mat *vir, struct sparse_mat *occ_trans, struct sparse_mat *vir_trans)
{
    struct sparse_vec work, work2;
    malloc_vec(fock->ncol, &work); malloc_vec(fock->ncol, &work2);
    // main solver loop
    int iter = 0;
    double metric;
    do
    {
int nnz = 0;
for(int i=0 ; i<occ->ncol ; i++)
{ nnz += occ->col[i].num; }
for(int i=0 ; i<vir->ncol ; i++)
{ nnz += vir->col[i].num; }
printf("average LMO size = %e\n", (double)nnz/(double)fock->ncol);
printf("1st LMO size = %d\n", occ->col->num);
        // calculate the projected Fock matrix
        struct sparse_mat fock_lmo;
        project_energy(fock, occ, energy_occ);
        project_energy(fock, vir, energy_vir);
        project_fock(cut_fock, fock, occ, vir_trans, &fock_lmo);
        double total_energy = 0.0;
        for(int i=0 ; i<occ->ncol ; i++)
        { total_energy += energy_occ[i]; }
        // calculate the convergence metric
        metric = 0.0;
nnz = 0;
        for(int i=0 ; i<fock_lmo.ncol ; i++)
        for(int j=0 ; j<fock_lmo.col[i].num ; j++)
        { metric += pow(fock_lmo.col[i].val[j],2); nnz++; }
        metric = sqrt(metric)/(double)fock_lmo.ncol;
printf("nnz = %d\n",nnz);
        printf("total energy @ iter %d = %e (%e)\n", ++iter, total_energy, metric);
        // eliminate couplings from the projected Fock matrix
        for(int i=0 ; i<fock_lmo.ncol ; i++)
        for(int j=0 ; j<fock_lmo.col[i].num ; j++)
        {
            // calculate Jacobi rotation angle [https://en.wikipedia.org/wiki/Jacobi_rotation]
            int k = fock_lmo.col[i].row[j];
            double coupling = fock_lmo.col[i].val[j];
            double beta = 0.5*(energy_vir[k] - energy_occ[i])/coupling;
            double rot = beta/(beta*beta + fabs(beta)*sqrt(beta*beta + 1.0));
            jacobi(cut_lmo, rot, i, k, occ, vir, occ_trans, vir_trans, &work, &work2);
        }
        free_mat(&fock_lmo);
    } while (metric > tol);
    free_vec(&work); free_vec(&work2);
}

// stub main to read from mozyme.dump-formatted files & run sparse pseudo-diagonalization
int main(int argc, char **argv)
{
    // read mozyme.dump file
    if(argc < 2)
    { printf("SYNTAX: <executable> <input file (mozyme.dump format)> <lmo_cut (default=1e-4)> <fock_cut (default=1e-4)> <tol (default=1e-2)>\n"); exit(1); }
    double cut_lmo = LMO_CUTOFF, cut_fock = FOCK_CUTOFF, tol = FOCK_TOLERANCE;
    FILE *data = fopen(argv[1], "r");
    if(argc > 2) { sscanf(argv[2], "%lf", &cut_lmo); }
    if(argc > 3) { sscanf(argv[3], "%lf", &cut_fock); }
    if(argc > 4) { sscanf(argv[4], "%lf", &tol); }
    int num_occ, num_vir;
    fscanf(data, "%d %d", &num_occ, &num_vir);
    struct sparse_mat fock, occ, vir;
    fock.ncol = num_occ + num_vir;
    occ.ncol = num_occ;
    vir.ncol = num_vir;
    read_mat(&occ, data);
    read_mat(&vir, data);
    read_mat(&fock, data);
    fclose(data);

    // construct transposes
    renormalize(&occ);
    renormalize(&vir);
    struct sparse_mat occ_trans, vir_trans;
    transpose_mat(fock.ncol, &occ, &occ_trans);
    transpose_mat(fock.ncol, &vir, &vir_trans);

    // run pseudo-diagonalization
    double *energy_occ = (double*)malloc(sizeof(double)*occ.ncol);
    double *energy_vir = (double*)malloc(sizeof(double)*vir.ncol);
    clock_t time = clock();
    pseudo_diag(cut_lmo, cut_fock, tol, energy_occ, energy_vir, &fock, &occ, &vir, &occ_trans, &vir_trans);
    time = clock() - time;

    // report timing, sparsity, & accuracy information
    double total_energy = 0.0;
    for(int i=0 ; i<num_occ ; i++)
    { total_energy += energy_occ[i]; }
    printf("final total energy = %e\n", total_energy);
    int nnz = 0;
    for(int i=0 ; i<occ.ncol ; i++)
    { nnz += occ.col[i].num; }
    for(int i=0 ; i<vir.ncol ; i++)
    { nnz += vir.col[i].num; }
    printf("average LMO size = %e\n", (double)nnz/(double)fock.ncol);
    printf("runtime = %e s\n", (double)time /(double)CLOCKS_PER_SEC);
    return 0;
}