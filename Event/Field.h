#ifndef __FIELD__
#define __FIELD__
#include <vector>
#include <string>
#include <utility>

class Field{
  public:
    Field(){ Init("/home/colin/GoogleDrive/437A/Code/potentialfield.txt"); };
    Field(std::string FieldFile){ Init(FieldFile); };

    std::vector<double> GetFieldAtPosition(double x, double y, double z);
    std::vector<double> GetFieldAtPosition(std::vector<double> _pos)
        { return GetFieldAtPosition(_pos[0],_pos[1],_pos[2]); };

    std::vector< std::vector<double> > GetField(){ return Potential; };
    std::vector<double> GetCoordinates(std::string _dir);

  protected:
    void Init(std::string FieldFile);

  private:
    std::vector< std::vector<double> > Potential;
    std::vector<double> CoordinatesX;
    std::vector<double> CoordinatesY;
    std::vector<double> CoordinatesZ;
};


#endif
