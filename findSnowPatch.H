    // 获取雪面的Patch名
    // get snow surface patch id
    std::vector<label> snowPatchList;
    const fvPatchList& patches = mesh.boundary();
    forAll(patches,i)
    {
        const std::string name = static_cast<std::string>(patches[i].name()); // 强制转换std::string
        std::smatch match;
        std::regex e("(.snow)");   // 匹配包含".snow"字符串的边界
        if(std::regex_search(name, match, e))
        {
            Info << "Recognized snow surface: \"" << patches[i].name() << "\"." << endl;
            snowPatchList.push_back(patches[i].index());
            phiWf.boundaryFieldRef()[patches[i].index()] = Zero; // 雪面处下落速度为0
        }
    }