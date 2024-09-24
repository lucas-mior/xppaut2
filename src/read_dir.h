#ifndef READ_DIR_H
#define READ_DIR_H
#include "integers.h"

typedef struct FileInfo {
    char **dirnames, **filenames;
    int32 nfiles;
    int32 ndirs;
} FileInfo;

#define MAXPATHLEN 512
extern char cur_dir[MAXPATHLEN];
extern FileInfo my_ff;

void free_finfo(FileInfo *ff);
int32 cmpstringp(const void *p1, const void *p2);
int32 get_fileinfo_tab(char *wild, char *direct, FileInfo *ff, char *wild2);
int32 get_fileinfo(char *wild, char *direct, FileInfo *ff);
int32 fil_count(char *direct, int32 *ndir, int32 *nfil, char *wild, int32 *mld,
                int32 *mlf);
int32 change_directory(char *path);
int32 get_directory(char *direct);
int32 is_directory(char *root, char *path);
void read_dir_make_full_path(char *root, char *filename, char *pathname);
int32 wild_match(char *string, char *pattern);

#endif
