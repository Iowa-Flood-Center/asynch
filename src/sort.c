#include "sort.h"

//Sorts a list of links by upstream area. Uses merge sort.
//Link** array: The list of links to be sorted.
//int size: The number of links in array.
void merge_sort(Link** array,unsigned int size,unsigned int idx)
{
	unsigned int increment,l,l_max,r,r_max,current,i;
	Link** tmp;
  
	increment = 1;
	tmp = (Link**) malloc(size * sizeof(Link*));
  
	while (increment < size)
	{
		l = 0;
		r = increment;
		l_max = r - 1;
		r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		current = 0;

		while (current < size)
		{
			while (l <= l_max && r <= r_max)
			{
				if (array[r]->params.ve[idx] > array[l]->params.ve[idx])
					tmp[current] = array[r++];
				else
					tmp[current] = array[l++];
				current++;
			}
	      
			while (r <= r_max)
				tmp[current++] = array[r++];
			while (l <= l_max)
				tmp[current++] = array[l++];

			l = r;
			r += increment;
			//l_max = r - 1;
			l_max = (r - 1 < size) ? r - 1: size - 1;
			r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		}
	    
		increment *= 2;
		for(i=0;i<size;i++)
			array[i] = tmp[i];
	}
	free(tmp);
}

//Sorts a list of links by distance. Uses merge sort.
//Link** array: The list of links to be sorted.
//int size: The number of links in array.
void merge_sort_distance(Link** sys,unsigned int* array,unsigned int size)
{
	unsigned int increment,l,l_max,r,r_max,current,i;
	unsigned int* tmp;
  
	increment = 1;
	tmp = (unsigned int*) malloc(size * sizeof(unsigned int));
  
	while (increment < size)
	{
		l = 0;
		r = increment;
		l_max = r - 1;
		r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		current = 0;

		while (current < size)
		{
			while (l <= l_max && r <= r_max)
			{
				if (sys[array[r]]->distance > sys[array[l]]->distance)
					tmp[current] = array[r++];
				else
					tmp[current] = array[l++];
				current++;
			}
	      
			while (r <= r_max)
				tmp[current++] = array[r++];
			while (l <= l_max)
				tmp[current++] = array[l++];

			l = r;
			r += increment;
			//l_max = r - 1;
			l_max = (r - 1 < size) ? r - 1: size - 1;
			r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		}
	    
		increment *= 2;
		for(i=0;i<size;i++)
			array[i] = tmp[i];
	}
	free(tmp);
}

//Sorts a two dimensional array by the ids entry (first column). Uses merge sort.
//int** array: The array to be sorted. The first column is assumed to hold the ids of each row. Must have two columns only.
//int size: The number of rows of array.
void merge_sort_ids(unsigned int** array,unsigned int size)
{
	unsigned int increment, l, l_max, r, r_max, current, i;
	unsigned int** tmp;

	increment = 1;
	tmp = (unsigned int**) malloc(sizeof(unsigned int*) * size);
	for(i=0;i<size;i++)	tmp[i] = malloc(sizeof(unsigned int) * 2);

	while (increment < size)
	{
		l = 0;
		r = increment;
		l_max = r - 1;
		r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		current = 0;

		while (current < size)
		{
			while (l <= l_max && r <= r_max)
			{
				if (array[r][0] < array[l][0])
				{
					tmp[current][0] = array[r][0];
					tmp[current][1] = array[r][1];
					r++;
				}
				else
				{
					tmp[current][0] = array[l][0];
					tmp[current][1] = array[l][1];
					l++;
				}
				current++;
			}

			while (r <= r_max)
			{
				tmp[current][0] = array[r][0];
				tmp[current][1] = array[r][1];
				current++;
				r++;
			}
			while (l <= l_max)
			{
				tmp[current][0] = array[l][0];
				tmp[current][1] = array[l][1];
				current++;
				l++;
			}

			l = r;
			r += increment;
			//l_max = r - 1;
			l_max = (r - 1 < size) ? r - 1: size - 1;
			r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		}

		increment *= 2;
		for (i = 0; i < size; i++)
		{
			array[i][0] = tmp[i][0];
			array[i][1] = tmp[i][1];
		}
	}
	for(i=0;i<size;i++)	free(tmp[i]);
	free(tmp);
}


//Sorts a one dimensional array by the entry. Uses merge sort.
//unsigned int* array: The array to be sorted.
//unsigned int size: The number of rows of array.
void merge_sort_1D(unsigned int* array,unsigned int size)
{
	unsigned int increment, l, l_max, r, r_max, current, i;
	unsigned int* tmp;

	increment = 1;
	tmp = (unsigned int*) malloc(sizeof(unsigned int*) * size);

	while (increment < size)
	{
		l = 0;
		r = increment;
		l_max = r - 1;
		r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		current = 0;

		while (current < size)
		{
			while (l <= l_max && r <= r_max)
			{
				if (array[r] < array[l])
				{
					tmp[current] = array[r];
					r++;
				}
				else
				{
					tmp[current] = array[l];
					l++;
				}
				current++;
			}

			while (r <= r_max)
			{
				tmp[current] = array[r];
				current++;
				r++;
			}
			while (l <= l_max)
			{
				tmp[current] = array[l];
				current++;
				l++;
			}

			l = r;
			r += increment;
			//l_max = r - 1;
			l_max = (r - 1 < size) ? r - 1: size - 1;
			r_max = (l_max + increment < size) ? l_max + increment : size - 1;
		}

		increment *= 2;
		for (i = 0; i < size; i++)
			array[i] = tmp[i];
	}

	free(tmp);
}

//Finds the location in sys of the link with id
unsigned int find_link_by_id(unsigned int id,Link** sys,unsigned int N)
{
	unsigned int k,max,min;

	k = N/2;
	max = N;
	min = 0;

	while(sys[k]->ID != id && max != min)
	{
		if(sys[k]->ID < id)	min = k+1;
		else			max = k;
		k = (max+min)/2;
	}

	if(sys[k]->ID != id)	return N+1;
	else			return k;
}

//Finds the location in sys of the link with id
unsigned int find_link_by_idtoloc(unsigned int id,unsigned int** id_to_loc,unsigned int N)
{
	unsigned int k,max,min;

	k = N/2;
	max = N;
	min = 0;

	while(max != min && id_to_loc[k][0] != id)
	{
		if(id_to_loc[k][0] < id)	min = k+1;
		else				max = k;
		k = (max+min)/2;
	}

	//if(id_to_loc[k][0] != id)	return N+1;
	if(max == min)			return N+1;
	else				return id_to_loc[k][1];
}

