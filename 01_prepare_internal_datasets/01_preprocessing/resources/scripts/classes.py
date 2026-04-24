"""
Classes for use with preprocessing pipeline scripts.
"""


class IndexKit:
    """
    Object class containing information about a sample index kit.
    """

    def __init__(self, file: str, lookup: dict, index_type: str) -> None:
        self.file = file
        self.lookup = lookup
        self.index_type = index_type
        assert isinstance(
            lookup, dict
        ), "Invalid attribute 'lookup'; must be an object of class 'dict'."
        assert index_type in {
            "dual",
            "single",
        }, f"Invalid attribute 'index_type': '{index_type}'; must be either 'dual' or 'single'."
        assert "index_name" in lookup.keys(), f"Invalid index kit CSV file: {file}"

    def __repr__(self) -> str:
        return f"IndexKit\n{self.index_type.capitalize()} index\nSource: {self.file}"

    def extract(self, name: str, index: int | None = None) -> str | tuple[str]:
        """
        Extract index sequence(s) based on index name.

        Arguments:
            ``name``: Index name.\n
            ``index``: Index sequence number (if more than 1 sequence associated with name).

        Returns:
            Index sequence (if ``index`` is specified) or tuple of index sequences
            (if ``index`` is ``None``).
        """
        return self.lookup[name][index] if index is not None else self.lookup[name]

    def match(self, names: set[str], strict: bool = True) -> bool:
        """
        Check if index kit contains index names.

        Arguments:
            ``names``: Set of index names.\n
            ``strict``: Boolean indicating if index name matching should be performed strictly.
            (see below).

        Returns:
            ``True`` if all (``strict = True``) or any (``strict = False``) index names in ``names``
            are present in index kit; otherwise ``False``.
        """
        names_match = (x in self.lookup.keys() for x in names)
        return all(names_match) if strict else any(names_match)
