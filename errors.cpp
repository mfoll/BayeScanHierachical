#include <fstream>
#include <stdlib.h>

//////////////////////////////////
// check file exist
/////////////////////////////////
inline void assure(std::ifstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}
inline void assure(std::ofstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}
inline void assure(std::fstream& in,
                   const char* filename = "")
{
    using namespace std;
    if (!in)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
}



