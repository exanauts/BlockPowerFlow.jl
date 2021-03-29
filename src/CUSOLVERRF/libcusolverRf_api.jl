# Julia wrapper for header: cusolverRf.h
# Automatically generated using Clang.jl

import CUDA: cuComplex, cuDoubleComplex
import CUDA.CUSOLVER: cusolverStatus_t
import CUDA.CUSPARSE: cusparseMatDescr_t

function cusolverRfCreate(handle)
    ccall((:cusolverRfCreate, CUSOLVER.libcusolver()), cusolverStatus_t, (Ptr{cusolverRfHandle_t},), handle)
end

function cusolverRfDestroy(handle)
    ccall((:cusolverRfDestroy, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t,), handle)
end

function cusolverRfGetMatrixFormat(handle, format, diag)
    ccall((:cusolverRfGetMatrixFormat, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{cusolverRfMatrixFormat_t}, Ptr{cusolverRfUnitDiagonal_t}), handle, format, diag)
end

function cusolverRfSetMatrixFormat(handle, format, diag)
    ccall((:cusolverRfSetMatrixFormat, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, cusolverRfMatrixFormat_t, cusolverRfUnitDiagonal_t), handle, format, diag)
end

function cusolverRfSetNumericProperties(handle, zero, boost)
    ccall((:cusolverRfSetNumericProperties, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Cdouble, Cdouble), handle, zero, boost)
end

function cusolverRfGetNumericProperties(handle, zero, boost)
    ccall((:cusolverRfGetNumericProperties, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cdouble}, Ptr{Cdouble}), handle, zero, boost)
end

function cusolverRfGetNumericBoostReport(handle, report)
    ccall((:cusolverRfGetNumericBoostReport, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{cusolverRfNumericBoostReport_t}), handle, report)
end

function cusolverRfSetAlgs(handle, factAlg, solveAlg)
    ccall((:cusolverRfSetAlgs, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, cusolverRfFactorization_t, cusolverRfTriangularSolve_t), handle, factAlg, solveAlg)
end

function cusolverRfGetAlgs(handle, factAlg, solveAlg)
    ccall((:cusolverRfGetAlgs, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{cusolverRfFactorization_t}, Ptr{cusolverRfTriangularSolve_t}), handle, factAlg, solveAlg)
end

function cusolverRfGetResetValuesFastMode(handle, fastMode)
    ccall((:cusolverRfGetResetValuesFastMode, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{cusolverRfResetValuesFastMode_t}), handle, fastMode)
end

function cusolverRfSetResetValuesFastMode(handle, fastMode)
    ccall((:cusolverRfSetResetValuesFastMode, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, cusolverRfResetValuesFastMode_t), handle, fastMode)
end

function cusolverRfSetupHost(n, nnzA, h_csrRowPtrA, h_csrColIndA, h_csrValA, nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU, h_P, h_Q, handle)
    ccall((:cusolverRfSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, cusolverRfHandle_t), n, nnzA, h_csrRowPtrA, h_csrColIndA, h_csrValA, nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU, h_P, h_Q, handle)
end

function cusolverRfSetupDevice(n, nnzA, csrRowPtrA, csrColIndA, csrValA, nnzL, csrRowPtrL, csrColIndL, csrValL, nnzU, csrRowPtrU, csrColIndU, csrValU, P, Q, handle)
    ccall((:cusolverRfSetupDevice, CUSOLVER.libcusolver()), cusolverStatus_t, (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, cusolverRfHandle_t), n, nnzA, csrRowPtrA, csrColIndA, csrValA, nnzL, csrRowPtrL, csrColIndL, csrValL, nnzU, csrRowPtrU, csrColIndU, csrValU, P, Q, handle)
end

function cusolverRfResetValues(n, nnzA, csrRowPtrA, csrColIndA, csrValA, P, Q, handle)
    ccall((:cusolverRfResetValues, CUSOLVER.libcusolver()), cusolverStatus_t, (Cint, Cint, CuPtr{Cint}, CuPtr{Cint}, CuPtr{Cdouble}, CuPtr{Cint}, CuPtr{Cint}, cusolverRfHandle_t), n, nnzA, csrRowPtrA, csrColIndA, csrValA, P, Q, handle)
end

function cusolverRfAnalyze(handle)
    ccall((:cusolverRfAnalyze, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t,), handle)
end

function cusolverRfRefactor(handle)
    ccall((:cusolverRfRefactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t,), handle)
end

function cusolverRfAccessBundledFactorsDevice(handle, nnzM, Mp, Mi, Mx)
    ccall((:cusolverRfAccessBundledFactorsDevice, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), handle, nnzM, Mp, Mi, Mx)
end

function cusolverRfExtractBundledFactorsHost(handle, h_nnzM, h_Mp, h_Mi, h_Mx)
    ccall((:cusolverRfExtractBundledFactorsHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), handle, h_nnzM, h_Mp, h_Mi, h_Mx)
end

function cusolverRfExtractSplitFactorsHost(handle, h_nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, h_nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU)
    ccall((:cusolverRfExtractSplitFactorsHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), handle, h_nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, h_nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU)
end

function cusolverRfSolve(handle, P, Q, nrhs, Temp, ldt, XF, ldxf)
    ccall((:cusolverRfSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, CuPtr{Cint}, CuPtr{Cint}, Cint, CuPtr{Cdouble}, Cint, CuPtr{Cdouble}, Cint), handle, P, Q, nrhs, Temp, ldt, XF, ldxf)
end

function cusolverRfBatchSetupHost(batchSize, n, nnzA, h_csrRowPtrA, h_csrColIndA, h_csrValA_array, nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU, h_P, h_Q, handle)
    ccall((:cusolverRfBatchSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, cusolverRfHandle_t), batchSize, n, nnzA, h_csrRowPtrA, h_csrColIndA, h_csrValA_array, nnzL, h_csrRowPtrL, h_csrColIndL, h_csrValL, nnzU, h_csrRowPtrU, h_csrColIndU, h_csrValU, h_P, h_Q, handle)
end

function cusolverRfBatchResetValues(batchSize, n, nnzA, csrRowPtrA, csrColIndA, csrValA_array, P, Q, handle)
    ccall((:cusolverRfBatchResetValues, CUSOLVER.libcusolver()), cusolverStatus_t, (Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Cint}, Ptr{Cint}, cusolverRfHandle_t), batchSize, n, nnzA, csrRowPtrA, csrColIndA, csrValA_array, P, Q, handle)
end

function cusolverRfBatchAnalyze(handle)
    ccall((:cusolverRfBatchAnalyze, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t,), handle)
end

function cusolverRfBatchRefactor(handle)
    ccall((:cusolverRfBatchRefactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t,), handle)
end

function cusolverRfBatchSolve(handle, P, Q, nrhs, Temp, ldt, XF_array, ldxf)
    ccall((:cusolverRfBatchSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{Cdouble}, Cint, Ptr{Ptr{Cdouble}}, Cint), handle, P, Q, nrhs, Temp, ldt, XF_array, ldxf)
end

function cusolverRfBatchZeroPivot(handle, position)
    ccall((:cusolverRfBatchZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverRfHandle_t, Ptr{Cint}), handle, position)
end

# Julia wrapper for header: cusolverSp_LOWLEVEL_PREVIEW.h
# Automatically generated using Clang.jl


function cusolverSpCreateCsrluInfoHost(info)
    ccall((:cusolverSpCreateCsrluInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (Ptr{csrluInfoHost_t},), info)
end

function cusolverSpDestroyCsrluInfoHost(info)
    ccall((:cusolverSpDestroyCsrluInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (csrluInfoHost_t,), info)
end

function cusolverSpXcsrluAnalysisHost(handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
    ccall((:cusolverSpXcsrluAnalysisHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t), handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
end

function cusolverSpScsrluBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpScsrluBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpDcsrluBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpDcsrluBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpCcsrluBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpCcsrluBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpZcsrluBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpZcsrluBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpScsrluFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
    ccall((:cusolverSpScsrluFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Cfloat, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
end

function cusolverSpDcsrluFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
    ccall((:cusolverSpDcsrluFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Cdouble, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
end

function cusolverSpCcsrluFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
    ccall((:cusolverSpCcsrluFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Cfloat, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
end

function cusolverSpZcsrluFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
    ccall((:cusolverSpZcsrluFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Cdouble, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pivot_threshold, pBuffer)
end

function cusolverSpScsrluZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpScsrluZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrluInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpDcsrluZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpDcsrluZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrluInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpCcsrluZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpCcsrluZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrluInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpZcsrluZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpZcsrluZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrluInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpScsrluSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpScsrluSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrluInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpDcsrluSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrluSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrluInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpCcsrluSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrluSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrluInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpZcsrluSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrluSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrluInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpXcsrluNnzHost(handle, nnzLRef, nnzURef, info)
    ccall((:cusolverSpXcsrluNnzHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t), handle, nnzLRef, nnzURef, info)
end

function cusolverSpScsrluExtractHost(handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
    ccall((:cusolverSpScsrluExtractHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cvoid}), handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
end

function cusolverSpDcsrluExtractHost(handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
    ccall((:cusolverSpDcsrluExtractHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cvoid}), handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
end

function cusolverSpCcsrluExtractHost(handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
    ccall((:cusolverSpCcsrluExtractHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cvoid}), handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
end

function cusolverSpZcsrluExtractHost(handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
    ccall((:cusolverSpZcsrluExtractHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrluInfoHost_t, Ptr{Cvoid}), handle, P, Q, descrL, csrValL, csrRowPtrL, csrColIndL, descrU, csrValU, csrRowPtrU, csrColIndU, info, pBuffer)
end

function cusolverSpCreateCsrqrInfoHost(info)
    ccall((:cusolverSpCreateCsrqrInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (Ptr{csrqrInfoHost_t},), info)
end

function cusolverSpDestroyCsrqrInfoHost(info)
    ccall((:cusolverSpDestroyCsrqrInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (csrqrInfoHost_t,), info)
end

function cusolverSpXcsrqrAnalysisHost(handle, m, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
    ccall((:cusolverSpXcsrqrAnalysisHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cint}, Ptr{Cint}, csrqrInfoHost_t), handle, m, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
end

function cusolverSpScsrqrBufferInfoHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpScsrqrBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrqrInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpDcsrqrBufferInfoHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpDcsrqrBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrqrInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpCcsrqrBufferInfoHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpCcsrqrBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrqrInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpZcsrqrBufferInfoHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpZcsrqrBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrqrInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpScsrqrSetupHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpScsrqrSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, Cfloat, csrqrInfoHost_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpDcsrqrSetupHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpDcsrqrSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cdouble, csrqrInfoHost_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpCcsrqrSetupHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpCcsrqrSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, cuComplex, csrqrInfoHost_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpZcsrqrSetupHost(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpZcsrqrSetupHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, cuDoubleComplex, csrqrInfoHost_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpScsrqrFactorHost(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpScsrqrFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpDcsrqrFactorHost(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrqrFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpCcsrqrFactorHost(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrqrFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpZcsrqrFactorHost(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrqrFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpScsrqrZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpScsrqrZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpDcsrqrZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpDcsrqrZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpCcsrqrZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpCcsrqrZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpZcsrqrZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpZcsrqrZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpScsrqrSolveHost(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpScsrqrSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpDcsrqrSolveHost(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrqrSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpCcsrqrSolveHost(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrqrSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpZcsrqrSolveHost(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrqrSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrqrInfoHost_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpXcsrqrAnalysis(handle, m, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
    ccall((:cusolverSpXcsrqrAnalysis, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cint}, Ptr{Cint}, csrqrInfo_t), handle, m, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
end

function cusolverSpScsrqrBufferInfo(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpScsrqrBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrqrInfo_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpDcsrqrBufferInfo(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpDcsrqrBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrqrInfo_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpCcsrqrBufferInfo(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpCcsrqrBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrqrInfo_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpZcsrqrBufferInfo(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpZcsrqrBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrqrInfo_t, Ptr{Cint}, Ptr{Cint}), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpScsrqrSetup(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpScsrqrSetup, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, Cfloat, csrqrInfo_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpDcsrqrSetup(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpDcsrqrSetup, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cdouble, csrqrInfo_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpCcsrqrSetup(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpCcsrqrSetup, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, cuComplex, csrqrInfo_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpZcsrqrSetup(handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
    ccall((:cusolverSpZcsrqrSetup, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, cuDoubleComplex, csrqrInfo_t), handle, m, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, mu, info)
end

function cusolverSpScsrqrFactor(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpScsrqrFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpDcsrqrFactor(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrqrFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpCcsrqrFactor(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrqrFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpZcsrqrFactor(handle, m, n, nnzA, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrqrFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, nnzA, b, x, info, pBuffer)
end

function cusolverSpScsrqrZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpScsrqrZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfo_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpDcsrqrZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpDcsrqrZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfo_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpCcsrqrZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpCcsrqrZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfo_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpZcsrqrZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpZcsrqrZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrqrInfo_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpScsrqrSolve(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpScsrqrSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpDcsrqrSolve(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrqrSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpCcsrqrSolve(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrqrSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpZcsrqrSolve(handle, m, n, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrqrSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrqrInfo_t, Ptr{Cvoid}), handle, m, n, b, x, info, pBuffer)
end

function cusolverSpCreateCsrcholInfoHost(info)
    ccall((:cusolverSpCreateCsrcholInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (Ptr{csrcholInfoHost_t},), info)
end

function cusolverSpDestroyCsrcholInfoHost(info)
    ccall((:cusolverSpDestroyCsrcholInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (csrcholInfoHost_t,), info)
end

function cusolverSpXcsrcholAnalysisHost(handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
    ccall((:cusolverSpXcsrcholAnalysisHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t), handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
end

function cusolverSpScsrcholBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpScsrcholBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpDcsrcholBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpDcsrcholBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpCcsrcholBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpCcsrcholBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpZcsrcholBufferInfoHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpZcsrcholBufferInfoHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpScsrcholFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpScsrcholFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpDcsrcholFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpDcsrcholFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpCcsrcholFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpCcsrcholFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpZcsrcholFactorHost(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpZcsrcholFactorHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpScsrcholZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpScsrcholZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpDcsrcholZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpDcsrcholZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpCcsrcholZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpCcsrcholZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfoHost_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpZcsrcholZeroPivotHost(handle, info, tol, position)
    ccall((:cusolverSpZcsrcholZeroPivotHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfoHost_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpScsrcholSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpScsrcholSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpDcsrcholSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrcholSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpCcsrcholSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrcholSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpZcsrcholSolveHost(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrcholSolveHost, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrcholInfoHost_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpCreateCsrcholInfo(info)
    ccall((:cusolverSpCreateCsrcholInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (Ptr{csrcholInfo_t},), info)
end

function cusolverSpDestroyCsrcholInfo(info)
    ccall((:cusolverSpDestroyCsrcholInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (csrcholInfo_t,), info)
end

function cusolverSpXcsrcholAnalysis(handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
    ccall((:cusolverSpXcsrcholAnalysis, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t), handle, n, nnzA, descrA, csrRowPtrA, csrColIndA, info)
end

function cusolverSpScsrcholBufferInfo(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpScsrcholBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpDcsrcholBufferInfo(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpDcsrcholBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpCcsrcholBufferInfo(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpCcsrcholBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpZcsrcholBufferInfo(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
    ccall((:cusolverSpZcsrcholBufferInfo, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cint}, Ptr{Cint}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, internalDataInBytes, workspaceInBytes)
end

function cusolverSpScsrcholFactor(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpScsrcholFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpDcsrcholFactor(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpDcsrcholFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpCcsrcholFactor(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpCcsrcholFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpZcsrcholFactor(handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
    ccall((:cusolverSpZcsrcholFactor, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Cint, cusparseMatDescr_t, Ptr{cuDoubleComplex}, Ptr{Cint}, Ptr{Cint}, csrcholInfo_t, Ptr{Cvoid}), handle, n, nnzA, descrA, csrValA, csrRowPtrA, csrColIndA, info, pBuffer)
end

function cusolverSpScsrcholZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpScsrcholZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpDcsrcholZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpDcsrcholZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpCcsrcholZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpCcsrcholZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Cfloat, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpZcsrcholZeroPivot(handle, info, tol, position)
    ccall((:cusolverSpZcsrcholZeroPivot, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Cdouble, Ptr{Cint}), handle, info, tol, position)
end

function cusolverSpScsrcholSolve(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpScsrcholSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, csrcholInfo_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpDcsrcholSolve(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpDcsrcholSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, csrcholInfo_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpCcsrcholSolve(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpCcsrcholSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, csrcholInfo_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpZcsrcholSolve(handle, n, b, x, info, pBuffer)
    ccall((:cusolverSpZcsrcholSolve, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, csrcholInfo_t, Ptr{Cvoid}), handle, n, b, x, info, pBuffer)
end

function cusolverSpScsrcholDiag(handle, info, diag)
    ccall((:cusolverSpScsrcholDiag, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Ptr{Cfloat}), handle, info, diag)
end

function cusolverSpDcsrcholDiag(handle, info, diag)
    ccall((:cusolverSpDcsrcholDiag, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Ptr{Cdouble}), handle, info, diag)
end

function cusolverSpCcsrcholDiag(handle, info, diag)
    ccall((:cusolverSpCcsrcholDiag, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Ptr{Cfloat}), handle, info, diag)
end

function cusolverSpZcsrcholDiag(handle, info, diag)
    ccall((:cusolverSpZcsrcholDiag, CUSOLVER.libcusolver()), cusolverStatus_t, (cusolverSpHandle_t, csrcholInfo_t, Ptr{Cdouble}), handle, info, diag)
end
