Material property tag
---------------------------

**Material property tag** = <integer>

**Description/Usage**

This card is required for continuation problems when the specified Continuation Type is "MT" or "UM". It identifies the property of the relevant material to be used as the continuation parameter. <integer> is the 4-digit property tag number assigned to that property in the file "mm_mp_const.h" as follows:

::

    #define TAGC_<PROPERTY_NAME> xxxx

**Examples**

If the continuation parameter is heat capacity (TAGC_HEAT_CAPACITY = 1600), use:

::

    Material property tag = 1600

**Technical Discussion**

Although a card for the chosen material property is still required in the \*.mat file for the relevant material, the property value specified there is overwritten by the continuation parameter.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
