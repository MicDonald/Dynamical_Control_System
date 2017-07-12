#ifndef DATA3D_H_INCLUDED
#define DATA3D_H_INCLUDED

#include <array>
#include <vector>
#include <cmath>
#include <iterator>


using data3d = std::array<double, 3>;

auto Increment = [] (
		data3d& a,
		const data3d& b
	) noexcept
	{
		a[0] += b[0];
		a[1] += b[1];
		a[2] += b[2];
	};


auto Substract = [] (
		const data3d& to,
		const data3d& from
	) noexcept
	{
		data3d diff;
		diff[0] = to[0] - from[0];
		diff[1] = to[1] - from[1];
		diff[2] = to[2] - from[2];
		return diff;
	};


auto Dot = [] (
		const data3d& data1,
		const data3d& data2
	) noexcept
	{
		return data1[0]*data2[0] + data1[1]*data2[1] + data1[2]*data2[2];
	};


auto Square = [] (
		const data3d& data
	) noexcept
	{
		return Dot(data, data);
	};


auto Equivalent = [] (
		const data3d& data1,
		const data3d& data2
	) noexcept
	{
		return	fabs(data1[0] - data2[0]) < 1.e-5 &&
			fabs(data1[1] - data2[1]) < 1.e-5 &&
			fabs(data1[2] - data2[2]) < 1.e-5;
	};


auto MatchAny = [] (
		const data3d& data1,
		const std::vector<data3d>& dataArray
	) noexcept
	{
		for ( const auto& data2 : dataArray )
			if ( Equivalent(data1, data2) )
				return true;
		return false;
	};


auto MatchIndex = [] (
		const data3d& data1,
		const std::vector<data3d>& dataArray
	) noexcept
	{
		for (	auto iter = dataArray.begin();
			iter != dataArray.end();
			++iter
		)
			if ( Equivalent(data1, *iter) )
				return distance(dataArray.begin(), iter);
		return static_cast<decltype(distance(dataArray.begin(), dataArray.end()))>(
			dataArray.size()
		);
	};

#endif // DATA3D_H_INCLUDED
