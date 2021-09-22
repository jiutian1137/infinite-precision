#include <iostream>
#include <iomanip>
#include <clmagic/calculation/number/floating.h>
#include <clmagic/calculation/number/function.h>

int main() {
  using calculation::FloatX;
  using calculation::Float128;
  using calculation::Float256;
  using calculation::calculate_pi;

  using Number = double;
  //using Number = FloatX<52, 11>;// fmod have some problem

  /**
   * precision = pow(2,digit2) = pow(16, digit16)
   *             pow(2,digit2) = pow(pow(2,4),digit16)
   *                           = pow(2,digit16*4)
   * 
   * so precision = pow(16,digit2/4), first bit may be error.
  */
  double multiplier = 1.0/3;
  int digit = std::numeric_limits<double>::digits/4 - 2;//std::numeric_limits<double>::digits/4 - 1, correct for [A,E]
  Float256 A = trunc((calculate_pi<double>()*multiplier) * pow(16.0,digit)) * pow(16.0,-digit);
  Float256 B = trunc(calculate_pi<double>(digit,multiplier) * pow(16.0,digit)) * pow(16.0,-digit-digit);
  Float256 C = trunc(calculate_pi<Number>(digit+digit,multiplier) * pow(16.0,digit)) * pow(16.0,-digit-digit-digit);
  std::cout << "         " << A << std::endl;
  std::cout << "         " << B << std::endl;
  std::cout << "         " << C << std::endl;
  std::cout << ".......: " << A + B + C << std::endl;
  std::cout << "correct: " << calculate_pi<Float256>() * multiplier << std::endl;
  std::cout << "bits   : ";
  for (size_t i = 0; i <= log(pow(16.0,digit*3))/log(10); ++i)
    if (i == 0) 
      std::cout << "1.";
    else
      std::cout << "1";
  std::cout << std::endl;

  /*
         6.28318530717958623199592693708837032318115234375
         0.000000000000000244929359829470296308106771699014453293891113716540530731
         0.000000000000000000000000000000339137106414736320829868593878860267219992
.......: 6.28318530717958647692528676655900576839433877908528316248499257680775072
correct: 6.28318530717958647692528676655900576839433879875021164194988918461563280
bits   : 1.1111111111111111111111111111111111111111111
  */


  using calculation::Float512;

  std::cout << "\n..." << std::endl;
  Float512 D = trunc(calculate_pi<double>(digit+digit+digit,multiplier) * pow(16.0,digit))                               * pow(16.0,-digit-digit-digit-digit);
  Float512 E = trunc(calculate_pi<double>(digit+digit+digit+digit,multiplier) * pow(16.0,digit))                         * pow(16.0,-digit-digit-digit-digit-digit);
  Float512 F = trunc(calculate_pi<double>(digit+digit+digit+digit+digit,multiplier) * pow(16.0,digit))                   * pow(16.0,-digit-digit-digit-digit-digit-digit);
  Float512 G = trunc(calculate_pi<double>(digit+digit+digit+digit+digit+digit,multiplier) * pow(16.0,digit))             * pow(16.0,-digit-digit-digit-digit-digit-digit-digit);
  Float512 H = trunc(calculate_pi<double>(digit+digit+digit+digit+digit+digit+digit,multiplier) * pow(16.0,digit))       * pow(16.0,-digit-digit-digit-digit-digit-digit-digit-digit);
  Float512 I = trunc(calculate_pi<double>(digit+digit+digit+digit+digit+digit+digit+digit,multiplier) * pow(16.0,digit)) * pow(16.0,-digit-digit-digit-digit-digit-digit-digit-digit-digit);
  Float512 J = trunc(calculate_pi<double>(digit+digit+digit+digit+digit+digit+digit+digit+digit,multiplier) * pow(16.0,digit)) * pow(16.0,-digit-digit-digit-digit-digit-digit-digit-digit-digit-digit);
  std::cout << "         " << Float512(A) << std::endl;
  std::cout << "         " << Float512(B) << std::endl;
  std::cout << "         " << Float512(C) << std::endl;
  std::cout << "         " << D << std::endl;
  std::cout << "         " << E << std::endl;
  std::cout << "         " << F << std::endl;
  std::cout << "         " << G << std::endl;
  std::cout << "         " << H << std::endl;
  std::cout << "         " << I << std::endl;
  std::cout << "         " << J << std::endl;
  std::cout << ".......: " << Float512(A) + Float512(B) + Float512(C) + D + E + F + G + H + I + J << std::endl;
  std::cout << "correct: " << calculate_pi<Float512>() * multiplier << std::endl;
  std::cout << "bits   : ";
  for (size_t i = 0; i <= log(pow(16.0,digit*10))/log(10); ++i)
    if (i == 0)
      std::cout << "1.";
    else
      std::cout << "1";
  std::cout << std::endl;

  /*
         1.04719755119657520481268875300884246826171875
         0.000000000000022483210384072354283048265804795425026885169472734560258687
         0.000000000000000000000000002252843786884131411276004392354916112033608144
.......: 1.04719755119659768802307282761596930341165495670103127752438884659386683
correct: 1.04719755119659768802307282761596930341170805657549344054511920252526874
bits   : 1.111111111111111111111111111111111111111

...
         1.04719755119657520481268875300884246826171875
         0.000000000000022483210384072354283048265804795425026885169472734560258686542510986328125
         0.000000000000000000000000002252843786884131411276004392354916112033608144412016007765423215337963203097615405567921698093414306640625
         0.000000000000000000000000000000000000000053099874462155635587032449260969606538563536543610988395699408011628348681957019649528312001812713205595173
         0.000000000000000000000000000000000000000000000000000007385143323481913011388384407321194254138613294888218980109530555561206205708357805012009055526
         0.000000000000000000000000000000000000000000000000000000000000000000227938832880181841624656054635603695843607481740317992942132498874079687512722944
         0.000000000000000000000000000000000000000000000000000000000000000000000000000000002712112029833846199886313040443793665221739699662119698872457447959
         0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001213976844661655097770922538994799131260870291298219515
         0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000071052031154379618168069649830666549756851
         0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000119450901389389223663554762
.......: 1.04719755119659768802307282761596930341170805657549344054511920252526875078233014950502276635466797782066438810147027080768429092427848646669635273
correct: 1.04719755119659768802307282761596930341170805657549344054511920252526875078233014950502276635466797782066438810147027080768429092427858669896112848
bits   : 1.111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
  */

  return 0;
}