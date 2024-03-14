#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int extraMemoryAllocated;

// void *Alloc(size_t sz)
// {
// 	extraMemoryAllocated += sz;
// 	size_t* ret = Alloc(sizeof(size_t) + sz);
// 	*ret = sz;
// 	printf("Extra memory allocated, size: %ld\n", sz);
// 	return &ret[1];
// }

void *Alloc(size_t sz)
{
    // Allocate memory for the size of the requested block plus the size of size_t
    size_t* ret = (size_t*)malloc(sz + sizeof(size_t));
    if (ret == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    // Store the size of the allocated block at the beginning of the memory block
    *ret = sz;

    // Increment the global variable to keep track of the extra memory allocated
    extraMemoryAllocated += sz + sizeof(size_t);

    printf("Extra memory allocated, size: %ld\n", sz);
    
    // Return a pointer to the memory block right after the size information
    return (void*)(ret + 1);
}


void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}

// Function to swap the position of two elements
 
void swap(int* a, int* b)
{
 
    int temp = *a;
    *a = *b;
    *b = temp;
}
 
// To heapify a subtree rooted with node i
// which is an index in arr[].
// n is size of heap

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapify(int arr[], int N, int i)
{
    // Find largest among root,
    // left child and right child
 
    // Initialize largest as root
    int largest = i;
 
    // left = 2*i + 1
    int left = 2 * i + 1;
 
    // right = 2*i + 2
    int right = 2 * i + 2;
 
    // If left child is larger than root
    if (left < N && arr[left] > arr[largest])
 
        largest = left;
 
    // If right child is larger than largest
    // so far
    if (right < N && arr[right] > arr[largest])
 
        largest = right;
 
    // Swap and continue heapifying
    // if root is not largest
    // If largest is not root
    if (largest != i) {
 
        swap(&arr[i], &arr[largest]);
 
        // Recursively heapify the affected
        // sub-tree
        heapify(arr, N, largest);
    }
}
 
// Main function to do heap sort
void heapSort(int arr[], int N)
{
 
    // Build max heap
    for (int i = N / 2 - 1; i >= 0; i--)
 
        heapify(arr, N, i);
 
    // Heap sort
    for (int i = N - 1; i >= 0; i--) {
 
        swap(&arr[0], &arr[i]);
 
        // Heapify root element
        // to get highest element at
        // root again
        heapify(arr, i, 0);
    }
}


// Merges two subarrays of arr[]. 
// First subarray is arr[l..m] 
// Second subarray is arr[m+1..r] 
void merge(int pData[], int l, int m, int r) 
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 = r - m; 
  
    // Create temp arrays 
    int L[n1], R[n2]; 
  
    // Copy data to temp arrays 
    // L[] and R[] 
    for (i = 0; i < n1; i++) 
        L[i] = pData[l + i]; 
    for (j = 0; j < n2; j++) 
        R[j] = pData[m + 1 + j]; 
  
    // Merge the temp arrays back 
    // into arr[l..r] 
    // Initial index of first subarray 
    i = 0; 
  
    // Initial index of second subarray 
    j = 0; 
  
    // Initial index of merged subarray 
    k = l; 
    while (i < n1 && j < n2) { 
        if (L[i] <= R[j]) { 
            pData[k] = L[i]; 
            i++; 
        } 
        else { 
            pData[k] = R[j]; 
            j++; 
        } 
        k++; 
    } 
  
    // Copy the remaining elements 
    // of L[], if there are any 
    while (i < n1) { 
        pData[k] = L[i]; 
        i++; 
        k++; 
    } 
  
    // Copy the remaining elements of 
    // R[], if there are any 
    while (j < n2) { 
        pData[k] = R[j]; 
        j++; 
        k++; 
    } 
} 

// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated
void mergeSort(int pData[], int l, int r)
{
    if(l < r){
        // set mid to --> ( left + right-1 ) / 2
        int m = l + (r - l) / 2; 
        
        // Sort first and second halves 
        mergeSort(pData, l, m); 
        mergeSort(pData, m + 1, r); 

        merge(pData, l, m, r);
    }

}
// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n)
{
	int i, item, j;
	for (i = 1; i < n; i++)
	{
		item = pData[i];
		/* Move elements of arr[0..i-1], that are
		greater than key, to one position ahead
		of their current position */
		for(j=i-1; j>=0; j--)
		{
			if(pData[j]>item)
				pData[j+1] = pData[j];
			else
				break;
		}
		pData[j+1] = item;
	}
}

// implement bubble sort
// extraMemoryAllocated counts bytes of extra memory allocated
void bubbleSort(int* pData, int n)
{
	//printf("\nUsing Bubble sort\n\n");
	int i, j,temp;
	for (i = 0; i < n-1; i++){
		//printf("Iteration# %d\n",i+1);
		for (j = 0; j < n-i-1; j++){
			if (pData[j] > pData[j+1]) { //then swap
				temp=pData[j];
				pData[j]=pData[j+1];
				pData[j+1]=temp;
			}
		}
	}
}

// implement selection sort
// extraMemoryAllocated counts bytes of extra memory allocated
void selectionSort(int* pData, int n)
{
	for (int i = 0; i < n - 1; i++) {
        
		int min_idx = i;

        for (int j = i + 1; j < n; j++) {
            if (pData[j] < pData[min_idx]) {
                min_idx = j;
            }
        }
        // Swap the found minimum element with the first element
        int temp = pData[i];
		pData[i] = pData[min_idx];
        pData[min_idx] = temp;

    }
}


// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	int i, n, *data;
	*ppData = NULL;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);

		// void *Alloc(size_t sz)
		*ppData = (int *)Alloc(sizeof(int) * dataSz);
		// Implement parse data block
		if (*ppData == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(-1);
		}
		for (i=0;i<dataSz;++i)
		{
			fscanf(inFile, "%d ",&n);
			data = *ppData + i;
			*data = n;
		}

		fclose(inFile);
	}
	
	return dataSz;
}

// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{

	int i, sz = (dataSz > 100 ? dataSz - 100 : 0);
	int firstHundred = (dataSz < 100 ? dataSz : 100);
	// int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<firstHundred;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;

	int i;
    double cpu_time_used;
	char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	
	for (i=0;i<1;++i)
	{
		int *pDataSrc, *pDataCopy;
		// Scan Number of Elements and Then Create Space for An Arry of Pointers to those elements
		int dataSz = parseData(fileNames[i], &pDataSrc);
		// IF NO elements continue
		if (dataSz <= 0)
			continue;
		// Allocate Space for a Copy of the Source Input to Be modified
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);

		// Number of elements in dataset
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.10lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.10lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.10lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.10lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

        printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.10lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
