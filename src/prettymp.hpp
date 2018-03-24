#include <iostream>
#include <iomanip>
#include "fixedmpz.hpp"
#include "mpmatrix.hpp"

namespace momentmp {
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

    class PrecPrint : virtual public PrettyMpMatrixBase {
      private:
        fmp_shift_t shift;
        void pretty_print(std::ostream &os) const {
            //Capture the initial flags of the output stream
            std::ios::fmtflags initialFlags(os.flags());

            size_t counter = 0;
            for (auto &mp : matrix) {
                counter++;
                os << std::setprecision(shift) << std::setw(1.5 * shift) << std::setfill(' ') 
                    << std::left << mp.to_mpf() << ' ';
                if (counter % matrix.getDimension() == 0) {
                    os << '\n';
                }
            }

            //Restore the initial flags of the output streams
            os.flags(initialFlags);
            os << std::endl;
        }
      public:
        PrecPrint(MpMatrix &matrix, fmp_shift_t shift) : shift(shift),
            PrettyMpMatrixBase(matrix) {}
    };

    class DebugPrint : public PrecPrint {
      public:
        DebugPrint(MpMatrix &matrix) : PrettyMpMatrixBase(matrix), PrecPrint(matrix, 5) {}
    };
}