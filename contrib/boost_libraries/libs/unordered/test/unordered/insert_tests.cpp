
// Copyright 2006-2010 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// clang-format off
#include "../helpers/prefix.hpp"
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "../helpers/postfix.hpp"
// clang-format on

#include "../helpers/test.hpp"
#include "../objects/test.hpp"
#include "../helpers/random_values.hpp"
#include "../helpers/tracker.hpp"
#include "../helpers/equivalent.hpp"
#include "../helpers/invariants.hpp"
#include "../helpers/input_iterator.hpp"
#include "../helpers/helpers.hpp"

#include <iostream>

namespace insert_tests {

test::seed_t initialize_seed(243432);

template <class X>
void unique_insert_tests1(X*, test::random_generator generator)
{
    test::check_instances check_;

    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;
    typedef test::ordered<X> ordered;

    std::cerr << "insert(value) tests for containers with unique keys.\n";

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);

    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {

        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        std::pair<iterator, bool> r1 = x.insert(*it);
        std::pair<BOOST_DEDUCED_TYPENAME ordered::iterator, bool> r2 =
            tracker.insert(*it);

        BOOST_TEST(r1.second == r2.second);
        BOOST_TEST(*r1.first == *r2.first);

        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    test::check_equivalent_keys(x);
}

template <class X>
void equivalent_insert_tests1(X*, test::random_generator generator)
{
    std::cerr << "insert(value) tests for containers with equivalent keys.\n";

    test::check_instances check_;

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);
    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {
        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        BOOST_DEDUCED_TYPENAME X::iterator r1 = x.insert(*it);
        BOOST_DEDUCED_TYPENAME test::ordered<X>::iterator r2 =
            tracker.insert(*it);

        BOOST_TEST(*r1 == *r2);

        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    test::check_equivalent_keys(x);
}

template <class X> void insert_tests2(X*, test::random_generator generator)
{
    typedef BOOST_DEDUCED_TYPENAME test::ordered<X> tracker_type;
    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;
    typedef BOOST_DEDUCED_TYPENAME X::const_iterator const_iterator;
    typedef BOOST_DEDUCED_TYPENAME tracker_type::iterator tracker_iterator;

    std::cerr << "insert(begin(), value) tests.\n";

    {
        test::check_instances check_;

        X x;
        tracker_type tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator r1 = x.insert(x.begin(), *it);
            tracker_iterator r2 = tracker.insert(tracker.begin(), *it);
            BOOST_TEST(*r1 == *r2);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <=
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert(end(), value) tests.\n";

    {
        test::check_instances check_;

        X x;
        X const& x_const = x;
        tracker_type tracker = test::create_ordered(x);

        test::random_values<X> v(100, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            const_iterator r1 = x.insert(x_const.end(), *it);
            tracker_iterator r2 = tracker.insert(tracker.end(), *it);
            BOOST_TEST(*r1 == *r2);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <=
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert(pos, value) tests.\n";

    {
        test::check_instances check_;

        X x;
        const_iterator pos = x.begin();
        tracker_type tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            pos = x.insert(pos, *it);
            tracker_iterator r2 = tracker.insert(tracker.begin(), *it);
            BOOST_TEST(*pos == *r2);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <=
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert single item range tests.\n";

    {
        test::check_instances check_;

        X x;
        tracker_type tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            x.insert(it, test::next(it));
            tracker.insert(*it);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <=
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert range tests.\n";

    {
        test::check_instances check_;

        X x;

        test::random_values<X> v(1000, generator);
        x.insert(v.begin(), v.end());

        test::check_container(x, v);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert range with rehash tests.\n";

    {
        test::check_instances check_;

        X x;

        test::random_values<X> v(1000, generator);

        x.insert(*v.begin());
        x.clear();

        x.insert(v.begin(), v.end());

        test::check_container(x, v);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert input iterator range tests.\n";

    {
        test::check_instances check_;

        X x;

        test::random_values<X> v(1000, generator);
        BOOST_DEDUCED_TYPENAME test::random_values<X>::const_iterator
            begin = v.begin(),
            end = v.end();
        x.insert(test::input_iterator(begin), test::input_iterator(end));
        test::check_container(x, v);

        test::check_equivalent_keys(x);
    }

    std::cerr << "insert copy iterator range tests.\n";

    {
        test::check_instances check_;

        X x;

        test::random_values<X> v(1000, generator);
        x.insert(test::copy_iterator(v.begin()), test::copy_iterator(v.end()));
        test::check_container(x, v);

        test::check_equivalent_keys(x);
    }

    std::cerr << "insert copy iterator range test 2.\n";

    {
        test::check_instances check_;

        X x;

        test::random_values<X> v1(500, generator);
        test::random_values<X> v2(500, generator);
        x.insert(
            test::copy_iterator(v1.begin()), test::copy_iterator(v1.end()));
        x.insert(
            test::copy_iterator(v2.begin()), test::copy_iterator(v2.end()));

        test::check_equivalent_keys(x);
    }

    std::cerr << "insert various ranges.\n";

    {
        for (int i = 0; i < 100; ++i) {
            X x;
            test::ordered<X> tracker = test::create_ordered(x);

            test::random_values<X> v(1000, generator);

            for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                     v.begin();
                 it != v.end();) {
                BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                    x.bucket_count();
                float b = x.max_load_factor();

                BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator next =
                    it;
                for (std::size_t j = test::random_value(20); j > 0; ++j) {
                    ++next;
                    if (next == v.end()) {
                        break;
                    }
                }

                x.insert(it, next);
                tracker.insert(it, next);
                it = next;

                tracker.compare(x); // Slow, but I can't see any other way.

                if (static_cast<double>(x.size()) <=
                    b * static_cast<double>(old_bucket_count))
                    BOOST_TEST(x.bucket_count() == old_bucket_count);
            }

            test::check_equivalent_keys(x);
        }
    }
}

template <class X>
void unique_emplace_tests1(X*, test::random_generator generator)
{
    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;
    typedef test::ordered<X> ordered;

    std::cerr << "emplace(value) tests for containers with unique keys.\n";

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);

    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {

        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        std::pair<iterator, bool> r1 = x.emplace(*it);
        std::pair<BOOST_DEDUCED_TYPENAME ordered::iterator, bool> r2 =
            tracker.insert(*it);

        BOOST_TEST(r1.second == r2.second);
        BOOST_TEST(*r1.first == *r2.first);

        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    tracker.compare(x);
    test::check_equivalent_keys(x);
}

template <class X>
void equivalent_emplace_tests1(X*, test::random_generator generator)
{
    std::cerr << "emplace(value) tests for containers with equivalent keys.\n";

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);
    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {
        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        BOOST_DEDUCED_TYPENAME X::iterator r1 = x.emplace(*it);
        BOOST_DEDUCED_TYPENAME test::ordered<X>::iterator r2 =
            tracker.insert(*it);

        BOOST_TEST(*r1 == *r2);

        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    tracker.compare(x);
    test::check_equivalent_keys(x);
}

template <class X> void move_emplace_tests(X*, test::random_generator generator)
{
    std::cerr
        << "emplace(move(value)) tests for containers with unique keys.\n";

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);

    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {

        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        typename X::value_type value = *it;
        x.emplace(boost::move(value));
        tracker.insert(*it);
        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    tracker.compare(x);
    test::check_equivalent_keys(x);
}

template <class X> void default_emplace_tests(X*, test::random_generator)
{
#if !BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x5100))
    std::cerr << "emplace() tests.\n";
    bool is_unique = test::has_unique_keys<X>::value;

    X x;

    x.emplace();
    BOOST_TEST(x.size() == 1);
    x.emplace();
    BOOST_TEST(x.size() == (is_unique ? 1u : 2u));
    x.emplace();
    BOOST_TEST(x.size() == (is_unique ? 1u : 3u));

    typename X::value_type y;
    BOOST_TEST(x.count(test::get_key<X>(y)) == (is_unique ? 1u : 3u));
    BOOST_TEST(*x.equal_range(test::get_key<X>(y)).first == y);

    x.emplace(y);
    BOOST_TEST(x.size() == (is_unique ? 1u : 4u));
    BOOST_TEST(x.count(test::get_key<X>(y)) == (is_unique ? 1u : 4u));
    BOOST_TEST(*x.equal_range(test::get_key<X>(y)).first == y);

    x.clear();
    BOOST_TEST(x.empty());
    x.emplace(y);
    BOOST_TEST(x.size() == 1);
    x.emplace(y);
    BOOST_TEST(x.size() == (is_unique ? 1u : 2u));

    BOOST_TEST(x.count(test::get_key<X>(y)) == (is_unique ? 1u : 2u));
    BOOST_TEST(*x.equal_range(test::get_key<X>(y)).first == y);
#endif
}

template <class X> void map_tests(X*, test::random_generator generator)
{
    std::cerr << "map tests.\n";

    X x;
    test::ordered<X> tracker = test::create_ordered(x);

    test::random_values<X> v(1000, generator);
    for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it = v.begin();
         it != v.end(); ++it) {
        BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count = x.bucket_count();
        float b = x.max_load_factor();

        x[it->first] = it->second;
        tracker[it->first] = it->second;

        tracker.compare_key(x, *it);

        if (static_cast<double>(x.size()) <=
            b * static_cast<double>(old_bucket_count))
            BOOST_TEST(x.bucket_count() == old_bucket_count);
    }

    tracker.compare(x);
    test::check_equivalent_keys(x);
}

template <class X> void map_tests2(X*, test::random_generator generator)
{
    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;
    std::cerr << "insert_or_assign\n";

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            std::pair<iterator, bool> r =
                x.insert_or_assign(it->first, it->second);
            BOOST_TEST(*r.first == *it);

            tracker[it->first] = it->second;
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert_or_assign(begin)\n";

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator r = x.insert_or_assign(x.begin(), it->first, it->second);
            BOOST_TEST(*r == *it);

            tracker[it->first] = it->second;
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert_or_assign(end)\n";

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator r = x.insert_or_assign(x.end(), it->first, it->second);
            BOOST_TEST(*r == *it);

            tracker[it->first] = it->second;
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }

    std::cerr << "insert_or_assign(last)\n";

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);
        iterator last = x.begin();

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator r = x.insert_or_assign(last, it->first, it->second);
            BOOST_TEST(*r == *it);

            tracker[it->first] = it->second;
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);

            last = r;
        }

        tracker.compare(x);
        test::check_equivalent_keys(x);
    }
}

template <class X> void try_emplace_tests(X*, test::random_generator generator)
{
    std::cerr << "try_emplace(key, value)\n";

    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator pos = x.find(it->first);
            bool found = pos != x.end();

            std::pair<typename X::iterator, bool> r =
                x.try_emplace(it->first, it->second);
            if (found) {
                BOOST_TEST(pos == r.first);
                BOOST_TEST(!r.second);
            } else {
                BOOST_TEST(r.second);
            }
            BOOST_TEST_EQ(r.first->first, it->first);
            BOOST_TEST_EQ(r.first->second, it->second);

            tracker.insert(*it);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        test::check_equivalent_keys(x);
    }

    std::cerr << "try_emplace(begin(), key, value)\n";

    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator pos = x.find(it->first);
            bool found = pos != x.end();

            typename X::iterator r =
                x.try_emplace(r.begin(), it->first, it->second);
            if (found) {
                BOOST_TEST(pos == r);
            }
            BOOST_TEST_EQ(r->first, it->first);
            BOOST_TEST_EQ(r->second, it->second);

            tracker.insert(*it);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        test::check_equivalent_keys(x);
    }

    std::cerr << "try_emplace(end(), key, value)\n";

    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator pos = x.find(it->first);
            bool found = pos != x.end();

            typename X::iterator r =
                x.try_emplace(r.end(), it->first, it->second);
            if (found) {
                BOOST_TEST(pos == r);
            }
            BOOST_TEST_EQ(r->first, it->first);
            BOOST_TEST_EQ(r->second, it->second);

            tracker.insert(*it);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        test::check_equivalent_keys(x);
    }

    std::cerr << "try_emplace(pos, key, value)\n";

    typedef BOOST_DEDUCED_TYPENAME X::iterator iterator;

    {
        test::check_instances check_;

        X x;
        test::ordered<X> tracker = test::create_ordered(x);

        test::random_values<X> v(1000, generator);
        for (BOOST_DEDUCED_TYPENAME test::random_values<X>::iterator it =
                 v.begin();
             it != v.end(); ++it) {
            BOOST_DEDUCED_TYPENAME X::size_type old_bucket_count =
                x.bucket_count();
            float b = x.max_load_factor();

            iterator pos = x.find(it->first);
            bool found = pos != x.end();

            typename X::iterator r = x.try_emplace(pos, it->first, it->second);
            if (found) {
                BOOST_TEST(pos == r);
            }
            BOOST_TEST_EQ(r->first, it->first);
            BOOST_TEST_EQ(r->second, it->second);

            tracker.insert(*it);
            tracker.compare_key(x, *it);

            if (static_cast<double>(x.size()) <
                b * static_cast<double>(old_bucket_count))
                BOOST_TEST(x.bucket_count() == old_bucket_count);
        }

        test::check_equivalent_keys(x);
    }
}

// Some tests for when the range's value type doesn't match the container's
// value type.

template <class X>
void map_insert_range_test1(X*, test::random_generator generator)
{
    std::cerr << "map_insert_range_test1\n";

    test::check_instances check_;

    typedef test::list<std::pair<BOOST_DEDUCED_TYPENAME X::key_type,
        BOOST_DEDUCED_TYPENAME X::mapped_type> >
        list;
    test::random_values<X> v(1000, generator);
    list l(v.begin(), v.end());

    X x;
    x.insert(l.begin(), l.end());

    test::check_equivalent_keys(x);
}

template <class X>
void map_insert_range_test2(X*, test::random_generator generator)
{
    std::cerr << "map_insert_range_test2\n";

    test::check_instances check_;

    typedef test::list<std::pair<BOOST_DEDUCED_TYPENAME X::key_type const,
        test::implicitly_convertible> >
        list;
    test::random_values<boost::unordered_map<BOOST_DEDUCED_TYPENAME X::key_type,
        test::implicitly_convertible> >
        v(1000, generator);
    list l(v.begin(), v.end());

    X x;
    x.insert(l.begin(), l.end());

    test::check_equivalent_keys(x);
}

boost::unordered_set<test::movable, test::hash, test::equal_to,
    std::allocator<test::movable> >* test_set_std_alloc;
boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, std::allocator<test::object> >* test_multimap_std_alloc;

boost::unordered_set<test::object, test::hash, test::equal_to,
    test::allocator1<test::object> >* test_set;
boost::unordered_multiset<test::movable, test::hash, test::equal_to,
    test::allocator2<test::movable> >* test_multiset;
boost::unordered_map<test::movable, test::movable, test::hash, test::equal_to,
    test::allocator2<test::movable> >* test_map;
boost::unordered_multimap<test::object, test::object, test::hash,
    test::equal_to, test::allocator1<test::object> >* test_multimap;

using test::default_generator;
using test::generate_collisions;
using test::limited_range;

UNORDERED_TEST(unique_insert_tests1,
    ((test_set_std_alloc)(test_set)(test_map))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(equivalent_insert_tests1,
    ((test_multimap_std_alloc)(test_multiset)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(
    insert_tests2, ((test_multimap_std_alloc)(test_set)(test_multiset)(
                       test_map)(test_multimap))(
                       (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(unique_emplace_tests1,
    ((test_set_std_alloc)(test_set)(test_map))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(equivalent_emplace_tests1,
    ((test_multimap_std_alloc)(test_multiset)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(move_emplace_tests,
    ((test_set_std_alloc)(test_multimap_std_alloc)(test_set)(test_map)(
        test_multiset)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(default_emplace_tests,
    ((test_set_std_alloc)(test_multimap_std_alloc)(test_set)(test_map)(
        test_multiset)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(map_tests,
    ((test_map))((default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(
    map_tests2, ((test_map))((default_generator)(generate_collisions)))

UNORDERED_TEST(map_insert_range_test1,
    ((test_multimap_std_alloc)(test_map)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

UNORDERED_TEST(map_insert_range_test2,
    ((test_multimap_std_alloc)(test_map)(test_multimap))(
                   (default_generator)(generate_collisions)(limited_range)))

#if !defined(BOOST_NO_CXX11_HDR_INITIALIZER_LIST)

struct initialize_from_two_ints
{
    int a, b;

    friend std::size_t hash_value(initialize_from_two_ints const& x)
    {
        return static_cast<std::size_t>(x.a + x.b);
    }

    bool operator==(initialize_from_two_ints const& x) const
    {
        return a == x.a && b == x.b;
    }
};

UNORDERED_AUTO_TEST(insert_initializer_list_set)
{
    boost::unordered_set<int> set;
    set.insert({1, 2, 3, 1});
    BOOST_TEST_EQ(set.size(), 3u);
    BOOST_TEST(set.find(1) != set.end());
    BOOST_TEST(set.find(4) == set.end());

    boost::unordered_set<initialize_from_two_ints> set2;

#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 5))
    set2.insert({{1, 2}});
#else
    set2.insert({1, 2});
#endif
    BOOST_TEST(set2.size() == 1);
    BOOST_TEST(set2.find({1, 2}) != set2.end());
    BOOST_TEST(set2.find({2, 1}) == set2.end());

    set2.insert({{3, 4}, {5, 6}, {7, 8}});
    BOOST_TEST(set2.size() == 4);
    BOOST_TEST(set2.find({1, 2}) != set2.end());
    BOOST_TEST(set2.find({3, 4}) != set2.end());
    BOOST_TEST(set2.find({5, 6}) != set2.end());
    BOOST_TEST(set2.find({7, 8}) != set2.end());
    BOOST_TEST(set2.find({8, 7}) == set2.end());

    set2.insert({{2, 1}, {3, 4}});
    BOOST_TEST(set2.size() == 5);
    BOOST_TEST(set2.find({1, 2}) != set2.end());
    BOOST_TEST(set2.find({2, 1}) != set2.end());
    BOOST_TEST(set2.find({3, 4}) != set2.end());
    BOOST_TEST(set2.find({5, 6}) != set2.end());
    BOOST_TEST(set2.find({7, 8}) != set2.end());
    BOOST_TEST(set2.find({8, 7}) == set2.end());
}

#if !BOOST_WORKAROUND(BOOST_MSVC, == 1800)

UNORDERED_AUTO_TEST(insert_initializer_list_multiset)
{
    boost::unordered_multiset<std::string> multiset;
    // multiset.insert({});
    BOOST_TEST(multiset.empty());
    multiset.insert({"a"});
    BOOST_TEST_EQ(multiset.size(), 1u);
    BOOST_TEST(multiset.find("a") != multiset.end());
    BOOST_TEST(multiset.find("b") == multiset.end());
    multiset.insert({"a", "b"});
    BOOST_TEST(multiset.size() == 3);
    BOOST_TEST_EQ(multiset.count("a"), 2u);
    BOOST_TEST_EQ(multiset.count("b"), 1u);
    BOOST_TEST_EQ(multiset.count("c"), 0u);
}

#endif

UNORDERED_AUTO_TEST(insert_initializer_list_map)
{
    boost::unordered_map<std::string, std::string> map;
    // map.insert({});
    BOOST_TEST(map.empty());
    map.insert({{"a", "b"}, {"a", "b"}, {"d", ""}});
    BOOST_TEST_EQ(map.size(), 2u);
}

UNORDERED_AUTO_TEST(insert_initializer_list_multimap)
{
    boost::unordered_multimap<std::string, std::string> multimap;
    // multimap.insert({});
    BOOST_TEST(multimap.empty());
    multimap.insert({{"a", "b"}, {"a", "b"}, {"d", ""}});
    BOOST_TEST_EQ(multimap.size(), 3u);
    BOOST_TEST_EQ(multimap.count("a"), 2u);
}

#endif

struct overloaded_constructor
{
    overloaded_constructor(int x1_ = 1, int x2_ = 2, int x3_ = 3, int x4_ = 4)
        : x1(x1_), x2(x2_), x3(x3_), x4(x4_)
    {
    }

    int x1, x2, x3, x4;

    bool operator==(overloaded_constructor const& rhs) const
    {
        return x1 == rhs.x1 && x2 == rhs.x2 && x3 == rhs.x3 && x4 == rhs.x4;
    }

    friend std::size_t hash_value(overloaded_constructor const& x)
    {
        std::size_t hash = 0;
        boost::hash_combine(hash, x.x1);
        boost::hash_combine(hash, x.x2);
        boost::hash_combine(hash, x.x3);
        boost::hash_combine(hash, x.x4);
        return hash;
    }
};

UNORDERED_AUTO_TEST(map_emplace_test)
{
    boost::unordered_map<int, overloaded_constructor> x;

#if !BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x5100))
    x.emplace();
    BOOST_TEST(
        x.find(0) != x.end() && x.find(0)->second == overloaded_constructor());
#endif

    x.emplace(2, 3);
    BOOST_TEST(
        x.find(2) != x.end() && x.find(2)->second == overloaded_constructor(3));

    x.try_emplace(5);
    BOOST_TEST(
        x.find(5) != x.end() && x.find(5)->second == overloaded_constructor());
}

UNORDERED_AUTO_TEST(set_emplace_test)
{
    boost::unordered_set<overloaded_constructor> x;
    overloaded_constructor check;

#if !BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x5100))
    x.emplace();
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);
#endif

    x.clear();
    x.emplace(1);
    check = overloaded_constructor(1);
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);

    x.clear();
    x.emplace(2, 3);
    check = overloaded_constructor(2, 3);
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);

    x.clear();
    x.emplace(4, 5, 6);
    check = overloaded_constructor(4, 5, 6);
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);

    x.clear();
    x.emplace(7, 8, 9, 10);
    check = overloaded_constructor(7, 8, 9, 10);
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);
}

struct derived_from_piecewise_construct_t
    : boost::unordered::piecewise_construct_t
{
};

derived_from_piecewise_construct_t piecewise_rvalue()
{
    return derived_from_piecewise_construct_t();
}

struct convertible_to_piecewise
{
    operator boost::unordered::piecewise_construct_t() const
    {
        return boost::unordered::piecewise_construct;
    }
};

UNORDERED_AUTO_TEST(map_emplace_test2)
{
    boost::unordered_map<overloaded_constructor, overloaded_constructor> x;

    x.emplace(boost::unordered::piecewise_construct, boost::make_tuple(),
        boost::make_tuple());
    BOOST_TEST(
        x.find(overloaded_constructor()) != x.end() &&
        x.find(overloaded_constructor())->second == overloaded_constructor());

    x.emplace(
        convertible_to_piecewise(), boost::make_tuple(1), boost::make_tuple());
    BOOST_TEST(
        x.find(overloaded_constructor(1)) != x.end() &&
        x.find(overloaded_constructor(1))->second == overloaded_constructor());

    x.emplace(piecewise_rvalue(), boost::make_tuple(2, 3),
        boost::make_tuple(4, 5, 6));
    BOOST_TEST(x.find(overloaded_constructor(2, 3)) != x.end() &&
               x.find(overloaded_constructor(2, 3))->second ==
                   overloaded_constructor(4, 5, 6));

    derived_from_piecewise_construct_t d;
    x.emplace(d, boost::make_tuple(9, 3, 1), boost::make_tuple(10));
    BOOST_TEST(x.find(overloaded_constructor(9, 3, 1)) != x.end() &&
               x.find(overloaded_constructor(9, 3, 1))->second ==
                   overloaded_constructor(10));

    x.clear();

    x.try_emplace(overloaded_constructor());
    BOOST_TEST(
        x.find(overloaded_constructor()) != x.end() &&
        x.find(overloaded_constructor())->second == overloaded_constructor());

    x.try_emplace(1);
    BOOST_TEST(
        x.find(overloaded_constructor(1)) != x.end() &&
        x.find(overloaded_constructor(1))->second == overloaded_constructor());

    x.try_emplace(overloaded_constructor(2, 3), 4, 5, 6);
    BOOST_TEST(x.find(overloaded_constructor(2, 3)) != x.end() &&
               x.find(overloaded_constructor(2, 3))->second ==
                   overloaded_constructor(4, 5, 6));
}

UNORDERED_AUTO_TEST(set_emplace_test2)
{
    boost::unordered_set<
        std::pair<overloaded_constructor, overloaded_constructor> >
        x;
    std::pair<overloaded_constructor, overloaded_constructor> check;

    x.emplace(boost::unordered::piecewise_construct, boost::make_tuple(),
        boost::make_tuple());
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);

    x.clear();
    x.emplace(boost::unordered::piecewise_construct, boost::make_tuple(1),
        boost::make_tuple(2, 3));
    check =
        std::make_pair(overloaded_constructor(1), overloaded_constructor(2, 3));
    ;
    BOOST_TEST(x.find(check) != x.end() && *x.find(check) == check);
}
}

RUN_TESTS()
