#include <iostream>
#include <iomanip>
#include "mpmatrix.hpp"

namespace mpmatrix {
    class PrettyMpMatrixBase {
      protected:
        MpMatrix &matrix;
        virtual void pretty_print(std::ostream &os) const = 0;

      public:
        PrettyMpMatrixBase(MpMatrix &matrix) : matrix(matrix) {}
        friend std::ostream &operator<<(std::ostream &os, const PrettyMpMatrixBase &matrix) {
            matrix.pretty_print(os);
            return os;
        }
    };

    class SetPrecPrint : virtual public PrettyMpMatrixBase {
      private:
        mp_bitcnt_t scale;
        void pretty_print(std::ostream &os) const {
            //Capture the initial flags of the output stream
            std::ios::fmtflags initialFlags(os.flags());

            size_t counter = 0;
            for (auto &mp : matrix) {
                counter++;
                os << std::setfill(' ') << std::left << (mp() >> scale);
                if (counter % matrix.getDimension() == 0) {
                    os << '\n';
                }
            }

            //Restore the initial flags of the output streams
            os.flags(initialFlags);
            os << std::endl;
        }
      public:
        SetPrecPrint(MpMatrix &matrix, mp_bitcnt_t scale) : scale(scale),
            PrettyMpMatrixBase(matrix) {}
    };

    class DebugPrint : public SetPrecPrint {
      public:
        DebugPrint(MpMatrix &matrix) : PrettyMpMatrixBase(matrix), SetPrecPrint(matrix, 5) {}
    };
}