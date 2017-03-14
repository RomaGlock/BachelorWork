#include "stl_io.h"
namespace stl
{
  using namespace std;
  
  const Group& Group::operator= (const std::pair<std::string, Group>& p)
  {
    (*this)=p.second;
    return *this;
  }
  
  Group::Group(const std::pair<std::string, Group>& p)
  {
    (*this)=p.second;
  }
}
