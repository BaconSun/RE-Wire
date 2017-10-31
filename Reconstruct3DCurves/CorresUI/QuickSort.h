#pragma once

#include "LibParty.h"

class Sorting {
public:
	static void QuickSort(vector<float> sum, vector<int> &vec, int left, int right)
	{
		if (left == right)
		{
			return;
		}

		int index = Partition(sum, vec, left, right);
		if (index > left)
		{
			QuickSort(sum, vec, left, index);
		}
		if (index < right)
		{
			QuickSort(sum, vec, index + 1, right);
		}
	}


	static int Partition(vector<float> sum, vector<int> &vec, int left, int right)
	{
		int i = left, j = right, x = vec[left];

		while (i < j)
		{
			while (i < j && sum[vec[j]] >= sum[x])
			{
				j--;
			}
			if (i < j)
			{
				vec[i++] = vec[j];
			}

			while (i < j && sum[vec[i]] < sum[x])
			{
				i++;
			}
			if (i < j)
			{
				vec[j--] = vec[i];
			}
		}
		vec[i] = x;

		return i;
	}

	static void SelectPair(vector<Vector2f> dists, float connect_thre, vector<int> index, vector<int> &pairTag) {

	}
};
