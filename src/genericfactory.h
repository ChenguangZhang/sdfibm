#pragma once

#include <unordered_map>
#include <memory>
#include <iostream>

namespace sdfibm {

template <typename BaseType, typename ParaType>
class GenericFactory
{
public:
    using Creator = std::function<std::unique_ptr<BaseType>(const ParaType&)>;
    using CreatorMap = std::unordered_map<std::string, Creator>;

public:
    static bool add(const std::string& name, Creator fn)
    {
        auto it = creators_.find(name);
        if(it == creators_.end())
        {
            creators_[name] = fn;
            return true;
        }
        return false;
    }

    static std::unique_ptr<BaseType> create(const std::string& name, const ParaType& para)
    {
        auto it = creators_.find(name);
        if(it != creators_.end())
        {
            return it->second(para);
        }
        else
        {
            throw std::runtime_error(std::string("Cannot create unrecognized object: " + name + '\n'));
        }
    }

    static void report(std::ostream& os = std::cout)
    {
        os << "Objects the factory can \"produce\":\n";
        int i = 1;
        for(const auto& it : creators_)
        {
            os << '[' << i++ << "] " << it.first << std::endl;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const GenericFactory<BaseType, ParaType>& factory)
    {
        factory.report(os);
        return os;
    }

private:
    static CreatorMap creators_;
};;

#define MAKESPECIALFACTORY(type, basetype, paratype) \
    using type##Factory = GenericFactory<basetype, paratype>; \
    template<> \
    type##Factory::CreatorMap type##Factory::creators_ = {};
}