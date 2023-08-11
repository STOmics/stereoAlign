#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:11 PM
# @Author  : zhangchao
# @File    : exceptions.py
# @Email   : zhangchao5@genomics.cn
class OptionalDependencyNotInstalled(ModuleNotFoundError):
    def __init__(self, exception, module_name=None):
        if module_name is None:
            module_name = exception.name
        self.message = (
            f"\n'{module_name}' is an optional dependency and not installed by default. "
            f"Please make sure you install it manually."
        )
        super().__init__(self.message)


class RLibraryNotFound(ModuleNotFoundError):
    def __init__(self, exception):
        self.message = f"\nproblem loading library: {exception}"
        super().__init__(self.message)
