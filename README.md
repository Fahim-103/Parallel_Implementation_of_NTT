# Parallel_Implementation_of_NTT

# Overview:
This repository contains the source codes for the parallel implementation of NTT on AI Engine. The procedure followed for the implementation is Daisuke Takahashi's "Six-Step NTT Algorithm".

* **Step-1:** In this step, we divide the 4096 elements in a 1024 X 4 matrix, meaning 1024 rows and 4 columns. We transpose the matrix resulting in a 4 X 1024 matrix, each row containing 1024 elements. For example, the first row would be (x0, x4, x8,.........,x4092).
* **Step-2:** Four 1024-point NTTs are executed in this step in four separate AI Engine processors. This step is again subdivided into a "Six-step NTT Algorithm" converting the 1024 point into a  4 X 256 matrix.
* **Step-3:** Twiddle multiplication is performed on each 1024-point NTT.
* **Step-4:** Matrix Transposition.
* **Step-5:** In this step, we perform 4-point NTTs for 1024 times.
* **Step-6:** Finally, we transpose the twiddle multiplied matrix and generate the output taking row-wise data.

Each of the four 1024-point NTT is computed in four distinct AI-Engine processors (kernels ntt_a, ntt_b, ntt_c, ntt_d ). The Output from each processor is then streamed to the final kernel processor (ntt_last), where we perform 4-point NTT for 1024 times. Each element of the four points comes from ntt_a, ntt_b, ntt_c, ntt_d kernels respectively.

## Calculation procedure of computing 1024-point NTT in each kernel.
For this computation, we follow NTT Radix-4 implementation using the "Six-Step NTT Algorithm". We divide the 1024-points into the following stages.

* **Stage0:** 1024 -> 4 X 265. We perform four-point NTT 256 times. After each NTT computation, we separate the first points in a 256-element chunk; the second points in another 256-element chunk; third and fourth points in two 256-point chunks. This way, we can perform the matrix Transposition stated in Step-4. After that, we multiply the elements with their corresponding pre-calculated twiddle values. All the chunks are stored in a single array (ntt_temp1), maintaining the separation.
 
* **Stage1:** 256 -> 4 X 64. We perform 4-point NTT 64 times on each 256-element chunk. And store the first, second, third, and fourth points in separate 64-element chunks. As the previous stage had four 256-element chunks, and every 256 chunks is distributed into four 64-element chunks hence, we can confirm we have executed operations on 4 X 4 X 64 = 1024 points. At the end of this stage, we perform corresponding twiddle multiplications.

* **Stage2:** 64 -> 4 X 16. Like the previous stage, we perform 4-point NTTs 16 times on each 64-element chunk. Next, we perform twiddle multiplications and store first, second, third, and fourth points accordingly.

* **Stage3+4:** 16 -> 4 X 4. We perform four 4-point NTTs on 16 elements and store the first, second, third, and fourth points in four separate chunks. Each chunk now consists of four elements. After twiddle multiplication, we perform 4-point NTT again on each four-element chunk. These two steps are covered in a single function in the source code.

As we are computing Radix-4 implementation of 1024-point NTT, we have five stages for generating 1024-point NTT (4^5 = 1024). In this document, the word "chunk" means "elements stored serially" in a temporary array.

![1024 Hierarchy](https://github.com/Fahim-103/Parallel_Implementation_of_NTT/assets/127331206/96c0f4dc-404e-4ec6-9fe5-f76e8937b5ad)


 
