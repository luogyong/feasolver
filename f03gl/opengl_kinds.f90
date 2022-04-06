MODULE opengl_kinds
USE, INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
PUBLIC

! Kind parameters
! Integer types:
INTEGER, PARAMETER :: GLbyte=C_SIGNED_CHAR, GLshort=C_SHORT,            &
    GLint=C_INT, GLsizei=C_INT, GLboolean=C_SIGNED_CHAR, GLenum=C_INT,  &
    GLbitfield=C_INT, GLcint=C_INT, GLubyte=C_SIGNED_CHAR, &
    GLushort=C_SHORT, GLuint=C_INT
! Real types:
INTEGER, PARAMETER :: GLdouble=C_DOUBLE, GLfloat=C_FLOAT, GLclampf=C_FLOAT, &
    GLclampd=C_DOUBLE

END MODULE opengl_kinds