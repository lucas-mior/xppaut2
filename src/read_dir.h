#ifndef READ_DIR_H
#define READ_DIR_H
#include "integers.h"

typedef struct {
    char **dirnames, **filenames;
    int32 nfiles;
    int32 ndirs;
} FILEINFO;

#define MAXPATHLEN 512
extern char cur_dir[MAXPATHLEN];
extern FILEINFO my_ff;

void free_finfo(FILEINFO *ff);
int32 cmpstringp(const void *p1, const void *p2);
int32 get_fileinfo_tab(char *wild, char *direct, FILEINFO *ff, char *wild2);
int32 get_fileinfo(char *wild, char *direct, FILEINFO *ff);
int32 fil_count(char *direct, int32 *ndir, int32 *nfil, char *wild, int32 *mld,
                int32 *mlf);
int32 change_directory(char *path);
int32 get_directory(char *direct);
int32 is_directory(char *root, char *path);
void MakeFullPath(char *root, char *filename, char *pathname);
int32 wild_match(char *string, char *pattern);

#endif
