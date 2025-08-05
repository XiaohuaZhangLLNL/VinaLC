#include <utility>
#include <sstream>
#include <string>
#include <iomanip>

#include <filesystem>
namespace fs = std::filesystem;

#ifndef VINALC_IO_UTILS
#define VINALC_IO_UTILS

namespace vinalc {
    namespace io_utils {
        /*! prepare output directory and perform
         *  necessary checks
         */
        void prep_outdir(const std::string &outdir,
                         const std::string &outfile);
    } //~ io_utils
} //~ vinalc
#endif
