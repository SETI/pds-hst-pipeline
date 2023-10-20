from typing import Any, Callable, Dict, List, TextIO, Tuple, TypeVar
import abc


class Lookup(abc.ABC):
    """
    A strategy for looking up a key value.  Subclasses may look in a
    single set of cards from a single file, or may look at multiple
    files, or do something completely different.
    """

    @abc.abstractmethod
    def __getitem__(self, key: str) -> Any:
        pass

    def get(self, key: str, default_value: Any) -> Any:
        try:
            return self.__getitem__(key)
        except KeyError:
            return default_value

    def keys(self, keys: List[str]) -> List[Any]:
        """Default implementation is to look up each key separated."""
        return list(map(self.__getitem__, keys))

    def dump_key(self, key: str, dump_file: TextIO) -> None:
        self.dump_keys([key], dump_file)

    @abc.abstractmethod
    def dump_keys(self, keys: List[str], dump_file: TextIO) -> None:
        pass

    @abc.abstractmethod
    def __repr__(self) -> str:
        pass


############################################################

CARD_SET = List[Dict[str, Any]]


class DictLookup(Lookup):
    """
    Look up a key value in a set of cards from one file.
    """

    def __init__(self, label: str, card_dicts: CARD_SET) -> None:
        self.label = label
        self.card_dicts = card_dicts

    def __getitem__(self, key: str) -> Any:
        return self.card_dicts[0][key]

    def dump_keys(self, keys: List[str], dump_file: TextIO) -> None:
        d = dict()
        for key in keys:
            try:
                val = self.__getitem__(key)
            except KeyError:
                val = None
            d[key] = val
        print(str(d), file=dump_file)

    def __repr__(self) -> str:
        return f"DictLookup({self.label})"


############################################################


class HduLookup(Lookup):
    """
    Look up a key value in a set of cards from one file.
    """

    def __init__(self, label: str, card_dict: Dict[str, Any]) -> None:
        self.label = label
        self.card_dict = card_dict

    def __getitem__(self, key: str) -> Any:
        return self.card_dict[key]

    def dump_keys(self, keys: List[str], dump_file: TextIO) -> None:
        d = dict()
        for key in keys:
            try:
                val = self.__getitem__(key)
            except KeyError:
                val = None
            d[key] = val
        print(str(d), file=dump_file)

    def __repr__(self) -> str:
        return f"HduLookup({self.label})"


def make_hdu_lookups(label: str, card_dicts: CARD_SET) -> List[Lookup]:
    return [
        HduLookup(f"{label}[{indx}]", card_dict)
        for indx, card_dict in enumerate(card_dicts)
    ]


# Merge Hdu0 and Hdu1. Some HST params data only exists in Hdu1 for some
# proposal id (ex: 13012), so when getting data HST params data, we will look
# at both Hdu0 and Hdu1.
def merge_two_hdu_lookups(hdu0: Lookup, hdu1: Lookup) -> Lookup:
    # Use data from Hdu0 as the primary, so if a param data exists in both
    # Hdu0 and Hdu1, we will use the one from hdu0.
    comb_dict = {**hdu1.card_dict, **hdu0.card_dict}  # type: ignore
    return HduLookup("Combine Hdu0 and Hdu1", comb_dict)


############################################################


class _DefaultDictLookup(Lookup):
    """
    Look up a key value in a set of cards from one file.
    """

    def __init__(self, label: str, card_dicts: CARD_SET) -> None:
        self.label = label
        self.card_dicts = card_dicts

    def __getitem__(self, key: str) -> Any:
        try:
            return self.card_dicts[0][key]
        except KeyError:
            return None

    def dump_keys(self, keys: List[str], f: TextIO) -> None:
        d = {key: self.__getitem__(key) for key in keys}
        print(f"{self.label}: {d}", file=f)

    def __repr__(self) -> str:
        return f"_DefaultDictLookup({self.label})"


############################################################

R = TypeVar("R")
D = TypeVar("D")


def _multi_lookup(f: Callable[[D], R], args: List[D]) -> R:
    """
    Runs a function over a list, returning the first successful result
    or raising the exception on the last try if there are no
    successes.
    """
    first_args = args[:-1]
    last_arg = args[-1]
    for arg in first_args:
        try:
            return f(arg)
        except KeyError:
            continue
    return f(last_arg)


################################


class MultiDictLookup(Lookup):
    def __init__(self, labeled_card_sets: List[Tuple[str, CARD_SET]]):
        self.labeled_card_sets = labeled_card_sets

    def __getitem__(self, key: str) -> Any:
        """Returns the first successful lookup."""

        def find_key(d: Tuple[str, CARD_SET]) -> Any:
            return d[1][0][key]

        return _multi_lookup(find_key, self.labeled_card_sets)

    def keys(self, keys: List[str]) -> List[Any]:
        """Returns the first successful lookup."""

        def find_keys(d: Tuple[str, CARD_SET]) -> List[Any]:
            return DictLookup(d[0], d[1]).keys(keys)

        return _multi_lookup(find_keys, self.labeled_card_sets)

    def dump_keys(self, keys: List[str], f: TextIO) -> None:
        for label, card_set in self.labeled_card_sets:
            _DefaultDictLookup(label, card_set).dump_keys(keys, f)

    def __repr__(self) -> str:
        return (
            f"MultiDictLookup({[label for label, card_set in self.labeled_card_sets]})"
        )
