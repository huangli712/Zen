!!! note

    This implementation of generic linked list in Fortran 90 was taken from

    ```text
    Jason R. Blevins's code
    journal: ACM Fortran Forum 28(3), 2-7, 2009.
    website: http://jblevins.org/research/generic-list.
    ```

    Of course, we have adapted the original code to fulfill our requirement.

This linked list is capable of storing data of an any type by using the generic programming techniques.

## Type

module

## Source

`src/m_linkedlist.f90`

## Usage

**(1)** Include linked list module support.

```fortran
use linkedlist
```

**(2)** Define user own data type.

```fortran
type data_t
     ! place your definition here
end type data_t
```

**(3)** Define pointer to data.

```fortran
type (data_t), pointer  :: data_ptr => null()
```

**(4)** Define pointer to list.

```fortran
type (list_t), pointer  :: list_ptr => null()
```

**(5)** Prepare data.

```fortran
allocate(data_ptr)
data_ptr%something = something
```

**(6)** Create a linked list.

```fortran
call list_init(list_ptr, transfer(data_ptr, list_d))
```

Here `list_d` is a public variable defined in linkedlist module.

!!! note

    To access the data stored in the nodes, we have to use the intrinsic `transfer()` subroutine.

**(7)** Insert new node.

```fortran
call list_insert(list_ptr, transfer(data_ptr, list_d))
```

`list_insert()` will always insert the new node after the given node, i.e., `list_ptr`. If you want to update the data for the given node, please use the following code:

```fortran
call list_put(list_ptr, transfer(data_ptr, list_d))
```

**(8)** Visit next node.

```fortran
curr => list_next(curr)
```

Here `curr` is a `list_t` type pointer.

**(9)** Retrieve data stored in the node.

```fortran
data_ptr  = transfer(list_get(curr), data_ptr)
something = data_ptr%something
```

Here `curr` is a `list_t` type pointer, it points to the current node.

**(10)** Free memory for this linked list.

```fortran
call list_free(list_ptr)
```
