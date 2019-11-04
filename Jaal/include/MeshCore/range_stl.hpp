#pragma once

#include <vector>
#include <set>
#include <algorithm>

namespace rstl
{
template<class T>
void sort(std::vector<T> &container)
{
    std::sort( container.begin(), container.end() );
}

template<class T>
void reverse(std::vector<T> &container)
{
    std::reverse( container.begin(), container.end() );
}

template<class T>
T min_element( const std::vector<T> &container)
{
    return *std::min_element(container.begin(), container.end() );
}

template<class T>
T max_element( const std::vector<T> &container)
{
    return *std::max_element(container.begin(), container.end() );
}

template<class T>
T mid_element( const std::vector<T> &container)
{
    size_t midindex = container.size()/2;
    return container[container.size()/2];
}

template<class T>
void remove_erase( std::vector<T> &v, const T &val)
{
    v.erase( std::remove( v.begin(), v.end(), val ), v.end() );
}

template<class T>
void remove( std::vector<T> &v, const T &val)
{
    std::remove(v.begin(), v.end(), val ), v.end();
}

template<class T>
void copy( const std::set<T> &src, std::vector<T> &dst)
{
    dst.clear();
    if( !src.empty() ) {
        dst.resize( src.size() );
        std::copy(src.begin(), src.end(), dst.begin() );
    }
}

template<class T>
bool contains( const std::vector<T> &src, const T &value)
{
    if( std::find( src.begin(), src.end(), value) == src.end() ) return false;
    return true;
}

template<class T>
void unique(std::vector<T> &data)
{
    auto last = std::unique(data.begin(), data.end());
    data.erase(last, data.end());
}

template<class T>
void set_difference(const std::set<T> &a, const std::set<T> &b, vector<T> &c)
{
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
}

template<class T>
void set_difference(const std::vector<T> &a, const std::vector<T> &b, vector<T> &c)
{
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
}

template<class T>
std::vector<T> set_difference(const std::vector<T> &a, const std::vector<T> &b)
{
    std::vector<T> c;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
    return c;
}

template<class T>
void set_intersection(const std::vector<T> &a, const std::vector<T> &b, vector<T> &c)
{
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
}

template<class T>
std::vector<T> set_intersection(const std::vector<T> &a, const std::vector<T> &b)
{
    std::vector<T> c;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
    return c;
}


}

